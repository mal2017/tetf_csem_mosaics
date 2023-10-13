library(AnnotationHub)
ah = AnnotationHub()

dm6_files <- query(ah, c("dm6","RepeatMasker"))


rpm <- dm6_files[["AH98983"]]

seqlevelsStyle(rpm) <- "NCBI"

rpm$name <- rpm$repName

rtracklayer::export(rpm, "ucsc_all_repeatmasker.ncbi_chrnames.bed")
