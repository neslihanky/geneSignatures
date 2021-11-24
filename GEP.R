pacman::p_load(MCPcounter,tidyverse,data.table)

expression <- fread("expressionmatrix.tsv") %>% 
  tibble::column_to_rownames("V1")

calculate_gep <-  function(countdf){
  houseKeep = c("ABCF1","NRDE2","G6PD",
                "OAZ1","POLR2A","SDHA",
                "STK11IP","TBC1D10B","TBP",
                "UBB","ZBTB34")
  GEPmarker = c("CCL5","CD27","CD274",
                "CD276","CD8A","CMKLR1",
                "CXCL9","CXCR6","HLA-DQA1",
                "HLA-DRB1","HLA-E","IDO1",
                "LAG3","NKG7","PDCD1LG2",
                "PSMB10","STAT1","TIGIT")
  # houseKeep <-  c("ENSG00000204574","ENSG00000119720","ENSG00000160211",
  #                 "ENSG00000104904","ENSG00000181222","ENSG00000073578",
  #                 "ENSG00000144589","ENSG00000169221","ENSG00000112592",
  #                 "ENSG00000170315","ENSG00000177125")
  # GEPmarker <-  c("ENSG00000161570","ENSG00000139193","ENSG00000120217",
  #                 "ENSG00000103855","ENSG00000153563","ENSG00000174600",
  #                 "ENSG00000138755","ENSG00000172215","ENSG00000196735",
  #                 "ENSG00000196126","ENSG00000204592","ENSG00000131203",
  #                 "ENSG00000089692","ENSG00000105374","ENSG00000197646",
  #                 "ENSG00000205220","ENSG00000115415","ENSG00000181847")
  GEPweight <-  c(0.008346,0.072293,0.042853,-0.0239,
                0.031021,0.151253,0.074135,0.004313,
                0.020091,0.058806,0.07175,0.060679,
                0.123895,0.075524,0.003734,0.032999,
                0.250229,0.084767)
  row.names(countdf) = gsub("\\..*$","",row.names(countdf)) # remove version from the Ensembl ID
  # follow the method from "IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade"
  hk.df <-  countdf[houseKeep,]
  hk.mean <-  colMeans(log10(hk.df+1),na.rm=T)
  gep.df <- sweep(log10(countdf[GEPmarker,]+1),2,hk.mean)
  outp <-  t(gep.df) %*% GEPweight
  colnames(outp) <-  "GEP"
  return(outp)
}

GEP <- calculate_gep(expression) %>% as.data.frame() %>%  tibble::rownames_to_column("Sample") 

write_tsv(GEP,"outputs/GEP_score.tsv")
