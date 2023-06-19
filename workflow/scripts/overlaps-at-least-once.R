library(tidyverse)
library(rtracklayer)
library(plyranges)

rpm <- import("~/work/tetf2/upstream/reference_insertions.bed")

subsetByOverlaps()


pan_tes <- read_tsv("~/work/tetf2/upstream/final-models.collected-info.tsv.gz") %>%
  filter(gene_symbol=="pan" & significant_x) %>%
  pull(feature.y)

chroi <- c("4","X","Y","2L","2R","3L","3R")
rpm.pan <- rpm %>% filter(name %in% pan_tes & seqnames %in%chroi)
rpm.other <- rpm %>% filter(!name %in% pan_tes & seqnames %in% chroi)

rpm.other %>% export("rpm.other.bed")
rpm.pan %>% export("rpm.pan.bed")
pan <- import("results/mosaics/w3l_pan_3/w3l_pan_3.mosaics.bed")


# if I include other chroms, worth changing this
all_ins_df <- rpm %>%
  filter(seqnames %in% c("4","X","Y","2L","2R","3L","3R")) %>%
  join_overlap_left(., pan) %>%
  as_tibble() %>%
  mutate(coex=name.x %in% pan_tes) %>%
  group_by(seqnames,start, end, name.x, coex) %>%
  summarise(bound = any(!is.na(name.y)))



all_ins_df %>%
  group_by(name.x, class=coex) %>%
  summarise(bound=any(bound),.groups = "drop") %>%
  group_by(class, bound) %>%
  tally() %>%
  mutate(bound=case_match(bound, T~"bound",
                          F~"unbound")) %>%
  mutate(class=case_match(class, T~"coex",
                          F~"ncoex")) %>%
  pivot_wider(names_from = bound, values_from = n, values_fill = 0) %>%
  dplyr::select(class,bound,unbound) %>%
  ungroup() %>%
  mutate(class=fct_relevel(class,c("coex","ncoex"))) %>%
  arrange(class) %>%
  column_to_rownames("class") %>%
  fisher.test()




pan %>%
  mutate(., pan.te = . %over% rpm.pan,
            other.te = . %over% rpm.other) %>%
  mutate(overlap = map2_chr(pan.te,other.te, .f =  function(x,y) {
    ifelse(x,"pan", ifelse(y, "other", "none"))
  })) %>%
  as_tibble() %>%
  ggplot(aes(overlap)) +
  geom_bar()
