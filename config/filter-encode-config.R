library(tidyverse)

st <- read_csv("config/sample_table_encode.csv")
sst <- read_csv("config/subsample_table_encode.csv")


filt_st <- st |> filter(target %in% c("pan","gro","arm","H3K9Me3","CG16779","Odj","CG7357","NfI","CG17802","vvl"))

filt_st <- bind_rows(filt_st, filter(st,sample_name %in% filt_st$input))

filt_sst <- filter(sst, sample_name %in% filt_st$sample_name)

write_csv(filt_st, "config/sample_table_encode.filt.csv")
write_csv(filt_sst, "config/subsample_table_encode.filt.csv")
