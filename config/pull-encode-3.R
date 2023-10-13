library(httr2)
library(jsonlite)
library(tidyverse)
library(furrr)
library(tictoc)

plan(multisession, workers = 4)

# ------------------------------
# get accessions I'm interested in
# ------------------------------
metadata_url <- "https://www.encodeproject.org/metadata/?status=released&replicates.library.biosample.donor.organism.scientific_name=Drosophila+melanogaster&type=Experiment"

metadata <- read_tsv(metadata_url)

# will work from fq files, only keeping first end to keep things consistent
# oddly, I found that some control ChIPs have a nonempty 'controlled by' section, which seems to have a totally random file in it
# e.g.ENCSR916PCV - ENCFF042HMR, which is apparently 'controlled by' ENCFF828CYL
# i replace these controlled by fields with NA to avoid issues later
metadata2 <- metadata |>
  filter(Assay %in% c("Control ChIP-seq", "TF ChIP-seq") & `File format` == "fastq" & (is.na(`Paired end`) | `Paired end` == 1)) |>
  mutate(`Controlled by` = if_else(Assay == "Control ChIP-seq", NA, `Controlled by`))

# remove low qual stuff
metadata3 <- metadata2 |> 
  filter(!`Audit NOT_COMPLIANT` %in% c("severe bottlenecking, poor library complexity,
                                      severe bottlenecking",
                                      "poor library complexity",
                                      "insufficient read depth",
                                      "nsufficient read depth, partially characterized antibody"))# ------------------------------
# adding requisite columns
# ------------------------------

# experiment to file lkup
experiment2file <- metadata3 |>
  dplyr::select(`Experiment accession`, `File accession`) |>
  distinct()

# prune some cols and add control file info
metadata4 <- metadata3 |>
  dplyr::select(`Experiment accession`,`Experiment target`,
                `File accession`, `Controlled by`, `Biological replicate(s)`,
                `Technical replicate(s)`) |>
  arrange(`Experiment accession`, `Biological replicate(s)`, `Technical replicate(s)`) |>
  mutate(`Controlled by` = str_extract(`Controlled by`,"(?<=files\\/).+?(?=\\/)")) |>
  mutate(target = str_extract(`Experiment target`,".+(?=-dmelanogaster)")) |>
  dplyr::select(-`Experiment target`) |>
  left_join(experiment2file, by = c(`Controlled by` =  "File accession")) |>
  dplyr::rename(`Experiment accession`="Experiment accession.x", Control.experiment = "Experiment accession.y")

# no technical replicates to worry about
stopifnot(nrow(metadata4 |> filter(str_detect(`Technical replicate(s)`,"_2")))==0)

# make a unique name that includes experiment and target and rep
metadata5 <- mutate(metadata4, sample = sprintf("%s_%s_rep%s",replace_na(target,"input"), `Experiment accession`, `Biological replicate(s)`)) |>
  dplyr::relocate(sample)

# sample name to experiment lookup (for adding control samples)
inputsample2experiment <- metadata5 |>
  filter(str_detect(sample,"^input_")) |>
  dplyr::select(sample, `Experiment accession`,`Biological replicate(s)`) |>
  dplyr::rename(input="sample") |>
  distinct()

# add input sample name that I will use
metadata6 <- left_join(metadata5, inputsample2experiment, by=c(Control.experiment="Experiment accession",
                                                  `Biological replicate(s)`="Biological replicate(s)"))

# encode annoyingly uses the same inputs for multiple biological replicates
# so here I fill in any missing input samples that were missed because I tried to join by biological replicate
# when making metadata6
metadata7 <-  split(filter(metadata6,!str_detect(sample,"input")), filter(metadata6,!str_detect(sample,"input"))$`Experiment accession`) |>
  map_df(~fill(.x, input, .direction="updown")) |>
  bind_rows(filter(metadata6,str_detect(sample,"input")))

# ------------
# Figured this out after running the full sample sheet...
# in cases where the biological replicates in input and IP don't match at all...
# we add in one of the replicates from the appropriate control experiment,
# even if the biological replicate doesn't match. This seems like the best option.
# only matters for ~200 samples/1600+
metadata7.5 <- filter(metadata7, !str_detect(sample,"input") & is.na(input)) |>
  #count(sample, sort = T) |> filter(n>1)
  left_join(inputsample2experiment, by=c(Control.experiment="Experiment accession"), relationship="many-to-many",suffix = c(".x",".y")) |>
  group_by(sample, `File accession`) |> # some samples have multiple files
  slice_head(n=1) |>
  ungroup() |>
  dplyr::select(-input.x,-`Biological replicate(s).y`) |>
  dplyr::rename(input="input.y", `Biological replicate(s)`="Biological replicate(s).x") |>
  bind_rows(filter(metadata7, !(!str_detect(sample,"input") & is.na(input))))
  
stopifnot(nrow(metadata7.5) == nrow(metadata7))

stopifnot(nrow(filter(metadata7.5, !str_detect(sample,"input") & is.na(input))) == 0)
# --------------

# uris in metadata file are obsoletes, plus not great sample info
# so I will pull this with the encode api
get_uri <- function(fac) {
  # This URL locates the ENCODE experiment
  message(fac)
  url = sprintf('https://www.encodeproject.org/files/%s/?frame=object', fac)
  
  req <- request(url) |> 
    req_headers("accept" = "application/json")
  
  resp <- req_perform(req)
  
  dat <- resp_body_json(resp)
  
  pluck(dat, "cloud_metadata","url")  
}


# get the public uri for every ENFF -  so i can dl with just curl/wget
metadata8 <- metadata7.5 |>
  mutate(accession = future_map_chr(`File accession`, get_uri))

# retrieve useful biosample info not included in metadata.tsv
get_biosample_info <- function(xac) {
  # This URL locates the ENCODE experiment
  message(xac)
  url = sprintf('https://www.encodeproject.org/experiments/%s/?frame=object', xac)
  
  req <- request(url) |> 
    req_headers("accept" = "application/json")
  
  resp <- req_perform(req)
  
  dat <- resp_body_json(resp)
  
  tibble(life_stage_age=pluck(dat, "life_stage_age"),
       perturbed=pluck(dat, "perturbed"),
       protein_tags = list(unique(pluck(dat, "protein_tags"))),
       simple_biosample_summary = pluck(dat, "simple_biosample_summary"))
}

# this is a REST GET for every experiment, but takes much less time than per file
metadata9 <- metadata8 |>
  dplyr::select(`Experiment accession`) |>
  distinct() |>
  mutate(biosample_info = future_map(`Experiment accession`, get_biosample_info)) |>
  left_join(metadata8, by = "Experiment accession")

# check to make sure theres only a single row in each biosample info so unnesting works
stopifnot(all(metadata9$biosample_info |> map_int(nrow))==1)

# unnest, add a column some of my pipelines use (SINGLE is enforced here because I am writing this to compare among
# other things mapq0 percentage, and I want a fairly level playing field across experiments)
metadata10 <- metadata9 |> 
  mutate(LIBRARY_LAYOUT="SINGLE") |>
  dplyr::rename(sample_name = "sample") |>
  dplyr::relocate(c(accession,LIBRARY_LAYOUT, input, biosample_info), .after = "sample_name") |>
  unnest(biosample_info) |>
  mutate(protein_tags =  map_chr(protein_tags, ~paste(map_chr(.x, ~pluck(.x,"name")),collapse = ";")))


# get histone chip stuff prior run
# originally made this for tetf_csem_mosaics, but I think I will use it here too
st3 <- read_csv("config/sample_table3.csv") |> dplyr::select(-`...5`)

st_histone <- st3 |> filter(str_detect(sample_name,"H3K|H4K") & str_starts(sample_name,"E"))

st_histone <- st_histone |> mutate(target = str_extract(sample_name,"H3K.+?(?=_)|H4K.+?(?=_)"))

st_histone <- bind_rows(st_histone, filter(st3, sample_name %in% st_histone$input))

full_st <- bind_rows(metadata10, st_histone)

# parens cause issues with shell scripts
full_st <- full_st |> mutate(sample_name = str_replace_all(sample_name,"[\\(\\)]","."))


full_st <- full_st |> mutate(sample_name =  str_replace_all(sample_name,"α","alpha"),
                  sample_name =  str_replace_all(sample_name,"β","beta"),
                  sample_name =  str_replace_all(sample_name,"γ","gamma"),
                  sample_name =  str_replace_all(sample_name,"δ","delta"))

                  

# now need to resolve "subsamples" -  some experiments have multiple fastqs per library
# it is a little unclear if these are different libraries from the same biological sample, resequenced libraries, or some mix of the two
full_st <- full_st |>
  group_by(sample_name) |>
  mutate(subsample_name = row_number()) |>
  ungroup() |>
  relocate(subsample_name, .after = sample_name)


sample_table <- full_st |>
  dplyr::select(sample_name, life_stage_age, simple_biosample_summary,target, input) |>
  distinct()

subsample_table <- full_st

write_csv(sample_table, "config/sample_table_encode.csv")
write_csv(subsample_table, "config/subsample_table_encode.csv")

# test
test_sample_table <- filter(sample_table, str_detect(sample_name,"^rn"))
test_sample_table <- bind_rows(test_sample_table, filter(sample_table, sample_name %in% test_sample_table$input))

test_subsample_table <- subsample_table |> filter(sample_name %in% test_sample_table$sample_name)

write_csv(test_sample_table, "config/test_sample_table_encode.csv")
write_csv(test_subsample_table, "config/test_subsample_table_encode.csv")
