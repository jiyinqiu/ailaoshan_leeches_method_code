library(tidyverse)
library(lulu)
library(knitr)
path_name <- file.path("/Volumes/JYQ/Ailaoshan/ailaoshan_leeches_method2019/analysis/16S_DAMe_SORT_outputs/Filter_min2PCRs_min9copies_16S/swarm_lulu/")
setwd(path_name)

assign(("matchlist_16S"), read.table("match_list_16S.txt", header = FALSE, as.is = TRUE, stringsAsFactors = FALSE))
assign(("16S_swarm_merge_otu"), read.table("16S_swarm_merge_otu.txt", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE))
assign("16S_swarm_merge_otu", get("16S_swarm_merge_otu") %>% column_to_rownames(var = "OTU")) 
assign("otutable_lulu", get("16S_swarm_merge_otu"))
assign("matchlist_lulu", get("matchlist_16S"))

curated_result <- lulu(otutable_lulu, matchlist_lulu)
curated_table <- rownames_to_column(curated_result$curated_table, var = "OTU")

write.table(curated_table, file = "16S_otu_table_swarm_lulu.txt", sep = "\t", row.names = FALSE, quote = FALSE)

rm(list=ls())
