pacman::p_load(data.table,tidyverse,fgsea)
setDTthreads(10)

## this script returns multiple gene signatures including Sia et al, inflammaotry, IFNG, exhausted t-cell by calculating GSVA scores using gene lists
expression <-fread("expressionmatrix.tsv") %>%   
  tibble::column_to_rownames("V1") 

Sia_immune_class <- fread("geneSets/Sia_et_al_immune_class_signature.txt") %>% 
  as.list()
names(Sia_immune_class) <- "ImmuneClass"
immune_signatures <- gmtPathways("geneSets/Immune_signatures.txt")

MVI_6_gene <- list(MVI_6_gene=c("ROS1", "UGT2B7", "FAS", "ANGPTL7", "GMNN", "MKI67"))

GS_full <- c(Sia_immune_class,immune_signatures,MVI_6_gene)


es.max<- GSVA::gsva(as.matrix(expression), GS_full, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)

write_tsv(es.max %>% t() %>% as.data.frame() %>% 
            tibble::rownames_to_column("sample") ,"outputs/GeneSignatures.tsv")

