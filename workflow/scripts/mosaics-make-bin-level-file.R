library(mosaics)

fraglen <- ifelse(exists("snakemake"),as.numeric(snakemake@params$fraglen),
               200)

binsize <- ifelse(exists("snakemake"),as.numeric(snakemake@params$binsize),
               50)

ibed <- ifelse(exists("snakemake"),snakemake@input$bed,
       "results/csem/embryo_0_8hr_pan_1.bed")

sizes <- ifelse(exists("snakemake"),snakemake@input$sizes,
               "~/BigData/dm6_peakseq_mappability/results/genome_filtered.sizes")

odir = ifelse(exists("snakemake"),
              snakemake@output$dir,
              "results/mosaics-bin-level/embryo_0_8hr_pan_1/")

dir.create(odir,recursive = T)

constructBins( infile=ibed,
  fileFormat="csem", 
  outfileLoc=odir,
  byChr=F, useChrfile=T, chrfile=sizes, excludeChr=NULL,
  PET=FALSE, fragLen=200, binSize=50, capping=0
)



  