library(tidyverse)
library(plyranges)
library(nullranges)
library(patchwork)
library(rtracklayer)
library(furrr)
library(tidyverse)
library(plyranges)
library(nullranges)
library(GenomicFeatures)
library(AnnotationDbi)
library(progress)

rpm <- import("~/work/tetf2/upstream/reference_insertions.bed") %>% sort()

peaks <- ifelse(exists("snakemake"),snakemake@input$bed, "results/mosaics-dedup/embryo_0_8hr_pan_1/embryo_0_8hr_pan_1.mosaics.bed") %>%
  import()

txdb <- "~/work/tetf2/results/resources/txdb" %>%
  loadDb()

# needed to compute the sequence lengths of the chromosomes
fa <- rtracklayer::import("~/work/tetf2/results/resources/genome.fasta")

L_s <- ifelse(exists("snakemake"),as.numeric(snakemake@params[["L_s"]]),1e5)
nseg <- ifelse(exists("snakemake"),as.numeric(snakemake@params[["nseg"]]),3)
R <- ifelse(exists("snakemake"),as.numeric(snakemake@params[["R"]]),100)
blockLength <- ifelse(exists("snakemake"),as.numeric(snakemake@params[["blockLength"]]),5e5)

# generate segmentation by gene density
g <- genes(txdb)

seqlengths(g) <- seqlengths(fa)[seqlevels(g)]
g <- keepStandardChromosomes(g, pruning.mode = "coarse")
g <- sortSeqlevels(g)
g <- sort(g)
names(g) <- NULL

seg_cbs <- segmentDensity(g, n = nseg, L_s = L_s, type = "cbs")

# prep the peaks for bootstrapping
seqlengths(peaks) <- seqlengths(fa)[seqlevels(peaks)]
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- sortSeqlevels(peaks)
peaks <- sort(peaks)
names(peaks) <- NULL
peaks$iter <- 0
peaks$name <- make.unique(peaks$name)

set.seed(2)
br <- bootRanges(peaks, blockLength, seg=seg_cbs, R=R,withinChrom=F)

rpm_gr_df <- split(rpm, rpm$name) %>% as.list() %>% enframe(name = "TE", value = "gr") 

# function to get our test statistic - the number of peaks in a set that overlap 
# a te
get_ol <- function(pks, TEs=rpm_gr_df) {
  mutate(TEs,n_overlaps = map_dbl(gr, ~{sum(pks %over% .x)})) %>%
    dplyr::select(-gr)
}

res <- c(peaks,br) %>% 
  split(.,.$iter) %>%
  as.list() %>%
  enframe(name = "boot", value = "gr") %>%
  mutate(stat = map(gr, get_ol))
  

res2 <- res %>%
  dplyr::select(-gr) %>%
  unnest(stat) %>%
  group_by(TE) %>%
  mutate(p = 1 - percent_rank(n_overlaps)) %>%
  ungroup()

write_tsv(res2, snakemake@output$tsv)



  # filter(boot == 0) %>%
  # ungroup() %>%
  # mutate(TE=fct_reorder(TE,p)) %>%
  # ggplot(aes(-log10(p),TE)) +
  # geom_point() +
  # geom_vline(xintercept = -log10(0.05),linetype="dashed",color="red") +
  # facet_wrap(~is.coex,scales = "free_y")







                                                     