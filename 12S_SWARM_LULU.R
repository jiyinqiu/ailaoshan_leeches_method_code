library(tidyverse)
library(lulu)
library(knitr)
path_name <- file.path("/Volumes/JYQ/Ailaoshan/ailaoshan_leeches_method2019/analysis/12S_DAMe_SORT_outputs/Filter_min2PCRs_min20copies_12S/swarm_lulu/")
setwd(path_name)

assign(("matchlist_12S"), read.table("match_list_12S.txt", header = FALSE, as.is = TRUE, stringsAsFactors = FALSE))
assign(("12S_swarm_merge_otu"), read.table("12S_swarm_merge_otu.txt", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE))
assign("12S_swarm_merge_otu", get("12S_swarm_merge_otu") %>% column_to_rownames(var = "OTU")) 
assign("otutable_lulu", get("12S_swarm_merge_otu"))
assign("matchlist_lulu", get("matchlist_12S"))

curated_result <- lulu(otutable_lulu, matchlist_lulu)
curated_table <- rownames_to_column(curated_result$curated_table, var = "OTU")

write.table(curated_table, file = "12S_otu_table_swarm_lulu.txt", sep = "\t", row.names = FALSE, quote = FALSE)

