
# Ok this is where it gets kind of weird, because csem won't normalize the wigs. Per the man page:
# >A wiggle plot representing the expected number of reads overlapping each position in the genome set can be generated from the sorted genome BAM file output. To generate the wiggle plot, run the 'csem-bam2wig' program on the '*.sorted.bam' files you have.
# See this for TMM based normalization https://www.biostars.org/p/413626/#414440
# I will do a simpler norm for now.

library(rtracklayer)

#ip_f <- "~/BigData/230315_tetf_chip_mosaics/results/csem_dedup/embryo_0_8hr_pan_1.bw"
ip_f <- snakemake@input[["ip"]]


#wce_f <- "~/BigData/230315_tetf_chip_mosaics/results/csem_dedup/embryo_0_8hr_pan_input_1.bw"
wce_f <- snakemake@input[["wce"]]

x <- ip_f |> import()

y <- wce_f |> import()

# these are binwidth of 1
#x |> width() |> hist()
#y |> width() |> hist()

x_tot <- sum(x$score)
y_tot <- sum(y$score)

x$score <- x$score / (x_tot/y_tot)

stopifnot(round(sum(x$score),1)==round(sum(y$score),1))

# alternatively, median scale
#med_x <- median(x$score)
#med_y <- median(y$score)
#x$score <- x$score/(med_x/med_y)
#stopifnot(median(x$score)==median(y$score))

export(x, snakemake@output$bw)