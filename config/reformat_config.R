# this reformats the config file from mal2017/tetf_csem_mosaics to work with nextflow
library(tidyverse)

x <- read_csv("config/sample_table3.csv")

x |>
    filter(!is.na(input)) |>
    filter(!str_detect(input,regex("input",ignore_case = T)))

x |>
    filter(input == sample_name)


x <- dplyr::select(x,sample=sample_name,control=input) |>
    mutate(fastq_1 = sprintf("results/basic_chip/fastq/%s.fastq.gz",sample),
        fastq_2="",
        antibody="") |>
    mutate(control = replace_na(control,"")) |>
    dplyr::select(sample,fastq_1,fastq_2,antibody, control)


# https://github.com/nf-core/chipseq/issues/313
# i recognize some samples may already have this, but the pipeline will fail if ANY don't....
# easier to add this and then strip the _REP1 from the names later
x <- x |> mutate(sample = paste0(sample,"_REP1"),
    control = if_else(control=="",control, paste0(control,"_REP1")))


x <- mutate(x, antibody = str_remove(str_remove(sample,regex("chip.+$|seq.+$|rep.+$|D\\.melano.+",ignore_case =T)),"_$|\\.$")) |> # nolint
    mutate(antibody =  if_else(control == "","",antibody)) |>
    arrange(desc(antibody))


# known to fail at peak calling
#x <- filter(x, !sample %in% c("H3K4me3_S2_cells_ChIP.seq_ChIP_Rep2_REP1","Su.Hw._16.24.seq_ChIP_Rep1_REP1"))


x_inputs <- filter(x,control=="")

x_min <- filter(x, str_detect(sample, "TCF|H3K|H4K") & control!="")

x <- filter(x_inputs, sample %in% x_min$control) |>
    bind_rows(x_min,y=_)

# needing if running a commit after tag 2.0.0
#x <- x |>
#    mutate(replicate = 1, control_replicate = 1) |>
#    dplyr::relocate(replicate, .before = "antibody") |>
#    dplyr::relocate(control_replicate, .after = "control")

x_test <- filter(x,str_detect(sample,"^TCF_0.8h.*Rep1"))




write_csv(x,"config/sample_table3.nfcore.csv")
write_csv(x_test,"config/sample_table3_test.nfcore.csv")

