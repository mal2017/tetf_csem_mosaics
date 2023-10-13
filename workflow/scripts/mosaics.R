library(mosaics)
library(tidyverse)

# filenames <- c(chip="results/mosaics-bin-level/embryo_0_8hr_pan_1/embryo_0_8hr_pan_1.bed_fragL200_bin50.txt",
#                input="results/mosaics-bin-level/embryo_0_8hr_pan_input_1/embryo_0_8hr_pan_input_1.bed_fragL200_bin50.txt",
#                M="~/BigData/dm6_peakseq_mappability/results/final-allchr/allchr_map_fragL200_bin50.txt",
#                GC="~/BigData/dm6_peakseq_mappability/results/final-allchr/allchr_GC_fragL200_bin50.txt",
#                N="~/BigData/dm6_peakseq_mappability/results/final-allchr/allchr_N_fragL200_bin50.txt")

# use params because parent rule outputs dir, i need to add to this string to get the file to use
filenames <- c(chip=snakemake@params$chip,
               input=snakemake@params$control,
               M=snakemake@input$M,
               GC=snakemake@input$GC,
               N=snakemake@input$N)

# ------------------------------------------------------------------------------
# get data + initial plotting
# ------------------------------------------------------------------------------
binTFBS <- readBins(type=names(filenames), fileName = filenames,dataType = "multi")

# ------------------------------------------------------------------------------
# fitting
# ------------------------------------------------------------------------------
fitTFBS <- mosaicsFit( binTFBS, bgEst="rMOM", analysisType="IO")

model2use <- names(which.min(c(`1S`=fitTFBS@bic1S,
            `2S`=fitTFBS@bic2S)))

# ------------------------------------------------------------------------------
# peak calling
# ------------------------------------------------------------------------------
peaks <- mosaicsPeak(fitTFBS, signalModel=model2use, FDR=0.1, maxgap=200, minsize=50, thres=10)

mosaics::export(peaks, type="bed", filename=snakemake@output$bed)

saveRDS(list(bins=binTFBS,fit=fitTFBS, peaks=peaks),snakemake@output$rds)

# peaks <- extractReads(peaks, chipPET=F,controlPET=F,
#              chipFile="results/csem/embryo_0_8hr_pan_1.bed",chipFileFormat="bed",chipFragLen=200,
#              controlFile="results/csem/embryo_0_8hr_pan_input_1.bed",controlFileFormat="bed", controlFragLen=200)
