#filter snps in biomarkers
disease_type_data  <- disease_type_daTa[,c(1,5,7,3)]


#disease_type_gene3 <- disease_type_gene[,c(1,5,7)]
annotated_snps4 = sqldf("
  SELECT *
  FROM annotated_snps3 d1 JOIN disease_type_data d2
  ON d1.gene = d2.gene
")
annotated_snps4 <- annotated_snps4[,-c(9)]
annotated_snps4 <-  annotated_snps4 %>% 
  mutate_all(funs(str_replace(., "no", ".")))
annotated_snps4<- annotated_snps4 %>% arrange(disease_type)
annotated_snps4 <- annotated_snps4[,c(11,2,1,9,10,3:8)]


