library(tidyverse)

# Goal of script is create a sample_table.csv with the following columns:
# sample_name,accession,LIBRARY_LAYOUT,input

df <- list(PRJNA63463 = "filereport_read_run_PRJNA63463_tsv.txt",
     PRJNA63467 = "filereport_read_run_PRJNA63467_tsv.txt") |>
  map_df(~set_tidy_names(read_tsv(.x),syntactic=T), .id="PROJECT")

# only take SE ChIP
#df$library_selection |> table()
df <- filter(df, library_selection == "ChIP" & library_layout == "SINGLE")

df <- dplyr::select(df, sample_name=sample_title, accession=fastq_ftp, LIBRARY_LAYOUT=library_layout) |>
  mutate(input="")

df <- df[!duplicated(df$sample_name),]

df$sample_name <- df$sample_name |>
  tidy_names(syntactic = T) |>
  str_replace("\\.+","\\.")


df <- df |>
  mutate(input = if_else(str_detect(sample_name,regex("input", ignore_case=T)),NA,input))


#search_terms <- c("CBP","nej","pan","TCF","gro","ctbp") |> 
#  paste(collapse="|")
#df |>
#  filter(str_detect(sample_name,regex(search_terms,ignore_case=T))) |>
#  print(n=Inf)


df |>
  distinct() %>% 
  write_csv("sample_table2_noinput.csv")
