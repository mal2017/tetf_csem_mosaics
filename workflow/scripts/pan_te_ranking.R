library(tidyverse)
library(plyranges)

# ----------------------------
# get TEs
message("importing te locations")
all_reps <- rtracklayer::import("/projectsp/f_cee53_1/ellison_lab/matt/tetf/subworkflows/tetf_refs/results/repeatmasker/genome.fasta.out.gff")
all_reps <- keepStandardChromosomes(all_reps, pruning.mode="coarse")
all_reps <- dropSeqlevels(all_reps,"Y", pruning.mode = "coarse")
all_reps$Target <- all_reps$Target |> str_extract('(?<=Motif:).+(?=")')

# this is super ubiquitous, so will exclude for speed
tes <- all_reps[!grepl("\\(", all_reps$Target) & !all_reps$Target == "INE-1",]

# ----------------------------
# get pan-assoc tes
message("importing lms")
pan_assoc <- read_tsv("/projectsp/f_cee53_1/ellison_lab/matt/tetf/subworkflows/tetf_lm_coex/results/linear_models/final-models.collected-info.tsv.gz") |> 
    filter(significant_x & gene_symbol == "pan") |> pull(feature.y)

tes <- tes |>
  mutate(grp = if_else(tes$Target %in% pan_assoc,"target","control")) |>
  mutate(id=paste(tes$Target,1:length(tes),sep="__"))

bw <- snakemake@input$bw

get_enrichment_per_insertion <- function(bw_fl, tes_gr=tes) {
  message(bw_fl)
  signal <- plyranges::read_bigwig(bw_fl)

  signal <- signal |>
    filter_by_overlaps(tes_gr) |>
    join_overlap_left(tes_gr,y=_) |>
    plyranges::group_by(id, grp, Target) |>
    reduce_ranges(enrichment=mean(score.y))

  return(as_tibble(signal))
  
}

message("running function")
res <- get_enrichment_per_insertion(bw, tes_gr = tes)

write_tsv(res, snakemake@output$tsv)