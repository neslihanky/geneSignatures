pacman::p_load(tidyverse,data.table)


expression <- fread("expressionmatrix.tsv") %>% 
  tibble::column_to_rownames("V1")


IMPRES <- fread("geneSets/IMPRES.txt") %>% 
  filter(ISORIGINAL=="GeneName") %>% select(-3) %>% 
  mutate_all(function(x) ifelse(x=="C10orf54","VSIR",x))## changed C10orf54 gene to VSIR 



all_impress_genes <- c(IMPRES$Gene1,IMPRES$Gene2) %>% unique


IMPRES_expression <- expression[all_impress_genes,]

IMPRES_df <- data.frame()
for (i in 1:15){
  newrow <- sapply(colnames(IMPRES_expression),function(x) ifelse(IMPRES_expression[IMPRES[i,1],x]<IMPRES_expression[IMPRES[i,2],x],1,0)) %>% 
    as.data.frame() %>% t
  IMPRES_df <- rbind(IMPRES_df,newrow)
}


IMPRES_df %>% colSums() %>% as.data.frame() %>% 
  tibble::rownames_to_column("Unique_ID") %>% rename(IMPRES=2) %>% 
  write_tsv("outputs/IMPRES_scores.txt")




