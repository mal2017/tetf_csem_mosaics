library(tidyverse)
library(rtracklayer)
library(plyranges)

pan_tes <- read_tsv("~/work/tetf2/upstream/final-models.collected-info.tsv.gz") %>%
  filter(significant_x) %>%
  filter(gene_symbol == "pan") %>%
  pull(feature.y) %>%
  unique()

tes <- import("~/work/tetf2/upstream/reference_insertions.bed")

peaks <- import("results/mosaics-dedup/embryo_0_8hr_pan_1/embryo_0_8hr_pan_1.mosaics.bed")

motifs <- import("~/work/tetf2/results/motifs/fimo_genome_wide/pan/fimo.gff")

s2r_leading_edge <- read_rds("~/work/tetf2/results/signatures/s2rplus_te_gsea.rds") %>%
  filter(RNAi=="pan") %>%
  pull(gsea) %>%
  .[[1]] %>%
  as_tibble() %>%
  pull(core_enrichment) %>%
  str_split("/") %>%
  .[[1]]

fls <- Sys.glob("results/enrichment/*tsv") %>%
  set_names(.,str_extract(.,"(?<=enrichment\\/).+(?=\\.te_enrich)"))

dat <- map_df(fls, read_tsv, .id="sample")

plot_enr <- function(x,nm="") {
  
  g <- x %>%
    filter(boot==0) %>%
    dplyr::select(-boot) %>%
    left_join(.,filter(dat,boot!=0), by=c("TE"),suffix=c(".original",".boots")) %>%
    mutate(TE=fct_reorder(TE,p.original)) %>%
    mutate(class = ifelse(TE %in% pan_tes,"coex.w.pan","other")) %>%
    filter(p.original < 0.05) %>%
    ggplot(aes(TE,n_overlaps.boots)) +
    geom_violin(scale = "area") +
    geom_point(aes(y=n_overlaps.original)) +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    facet_grid(~class,scales="free_x",space="free") +
    ggtitle(nm)

    return(g)
}

dat %>%
  nest(-sample) %>%
  deframe() %>%
  imap(plot_enr)


# ---------------------------
if using this in the future, be sure to exclude tes with fixed insertions
around each te, either here or upstream in the pipeline (230331)
hits <- dat %>%
  filter(sample == "embryo_0_8hr_pan_1" & boot == 0) %>%
  filter(p < 0.05) %>%
  dplyr::select(gs_name=sample, ensembl_gene = TE) %>%
  distinct()

library(clusterProfiler)

res_path <- "~/work/tetf2/results/deg/s2rplus.res.tsv.gz"
res <- read_tsv(res_path)

possibly_gsea <-  possibly(function(.x,.y) {
  set.seed(2022)
  GSEA(deframe(.x), TERM2GENE = hits,seed=2022,pvalueCutoff = 2, minGSSize = 5, eps = 0)
}, NULL)

gsea.tbl <- res %>%
  filter(comparison == "pan") %>%
  dplyr::select(comparison,feature,score=t) %>%
  arrange(-score) %>%
  nest(c(feature,score)) %>%
  #filter(comparison == "pan") %>%
  #mutate(gsea = future_map(data,possibly_gsea,.options=furrr_options(seed=2022))) %>%
  mutate(gsea = map2(data,comparison,possibly_gsea)) %>%
  mutate(gsea.tidy = map(gsea,as_tibble))

gsea.res <- gsea.tbl %>%
  dplyr::select(comparison,gsea,gsea.tidy) %>%
  mutate(kd = comparison) %>%
  unnest(gsea.tidy) %>%
  group_by(comparison) %>%
  mutate(pvRnk = dense_rank(-log10(pvalue)),
         nesRnk = dense_rank(abs(NES))) %>%
  ungroup() %>%
  relocate(pvRnk,nesRnk) %>%
  mutate(p.adjust = p.adjust(pvalue, method="BH")) %>%
  dplyr::select(RNAi=comparison,TE.set=ID,setSize,NES,pvalue, gsea)

enrichplot::gseaplot2(geneSetID = "embryo_0_8hr_pan_1", gsea.res$gsea[[1]],pvalue_table = T)
