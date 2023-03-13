###prepare data

#vcf_indel <- read.vcfR("demo_data/indel/control/c1_SRR12878708.vcf.bgz")
#head(vcf_indel)
#test <- readVcfFields('demo_data/indel/control/c1_SRR12878708.vcf.bgz', fields=c('CHROM','POS','REF','ALT'))
#part I : controls
#library("R.utils")
#df=gunzip("demo_data/indel/patients/p3_SRR12878719.vcf.bgz","demo_data/indel/patients/p3_SRR12878719.vcf")

#########patient

file_contenets_patients_indel <- list()
file_contenets_patients_indel <-file_paths_patients_indel %>%
  map(function(path){
    fread(path)
  })
#add each list item into the environment variable as a tibble.
# and set each list item to an environment variable
file_contenets_patients_indel%>% map(as_tibble) %>% 
  list2env(envir = .GlobalEnv)
# To join all files together into one data frame7
comined_patients_indel <-file_contenets_patients_indel %>% map(as_tibble) %>%
  reduce(full_join)
# delete NAs
# result : keep commun between all samples
comined_patients_indel <- na.omit(comined_patients_indel)
## identify snps
comined_patients_indel <- comined_patients_indel[,c(1:5)]
colnames(comined_patients_indel) <- c("chrom","POS_snp_patient","id_patient","REF_patient","ALT_patient")

###########"control

file_contenets_control_indel <- list()
file_contenets_control_indel <-file_paths_control_indel %>%
  map(function(path){
    fread(path)
  })
#add each list item into the environment variable as a tibble.
# and set each list item to an environment variable
file_contenets_control_indel%>% map(as_tibble) %>% 
  list2env(envir = .GlobalEnv)
# To join all files together into one data frame7
comined_control_indel <-file_contenets_control_indel %>% map(as_tibble) %>%
  reduce(full_join)
# delete NAs
# result : keep commun between all samples
comined_control_indel <- na.omit(comined_control_indel)
#
## identify snps
comined_control_indel <- comined_control_indel[,c(1:5)]
colnames(comined_control_indel) <- c("chrom","POS_snp_control","id_control","REF_control","ALT_control")

## add gene
##############################################################################################
chrom_snppos <- comined_control_indel[,c(1,2)]
chrom_snppos <-chrom_snppos %>%  mutate_all(funs(str_replace(., "chr", "")))
names(chrom_snppos)[2] <- "snp_position"
names(chrom_snppos)[1] <- "chrom"
#
gene_pos <- unique(annotation1[,c(5,6,7,9)] )#>>> filet genes
names(gene_pos)[1] <- "chrom"
names(gene_pos)[2] <- "start"
names(gene_pos)[3] <- "end"


annotated_indel = sqldf("
  SELECT *
  FROM gene_pos d1 JOIN chrom_snppos d2
  ON d1.chrom = d2.chrom
  AND d2.snp_position < d1.end
  AND d2.snp_position > d1.start
")
annotated_indel<-annotated_indel[,-c(5)]
#merge with snp info
comined_control_indel <-comined_control_indel %>%  mutate_all(funs(str_replace(., "chr", "")))
annotated_indel2 = sqldf("
  SELECT *
  FROM annotated_indel d1 JOIN comined_control_indel d2
  ON d1.chrom = d2.chrom
  AND d1.snp_position = d2.POS_snp_control
")
annotated_indel2 <- annotated_indel2[,-c(5,6)]
##############################################################################################""
chrom_snppos <- comined_patients_indel[,c(1,2)]
chrom_snppos <-chrom_snppos %>%  mutate_all(funs(str_replace(., "chr", "")))
names(chrom_snppos)[2] <- "snp_position"
names(chrom_snppos)[1] <- "chrom"
#
gene_pos <- unique(annotation1[,c(5,6,7,9)] )#>>> filet genes
names(gene_pos)[1] <- "chrom"
names(gene_pos)[2] <- "start"
names(gene_pos)[3] <- "end"


annotated_indel3 = sqldf("
  SELECT *
  FROM gene_pos d1 JOIN chrom_snppos d2
  ON d1.chrom = d2.chrom
  AND d2.snp_position < d1.end
  AND d2.snp_position > d1.start
")
annotated_indel3 <- annotated_indel3[,-c(5)]
#merge with snp info
comined_patients_indel <-comined_patients_indel %>%  mutate_all(funs(str_replace(., "chr", "")))
annotated_indel4 = sqldf("
  SELECT *
  FROM annotated_indel d1 JOIN comined_patients_indel d2
  ON d1.chrom = d2.chrom
  AND d1.snp_position = d2.POS_snp_patient
")
annotated_indel4 <- annotated_indel4[,-c(5,6)]
names(annotated_indel4)[4] <- "gene"
names(annotated_indel2)[4] <- "gene"
indels_patients <- annotated_indel4
indels_controls <-annotated_indel2
names(indels_patients)[5] <- "snp_position"
names(indels_controls)[5] <- "snp_position"
#### Part III: compare groups
compare_group_i <- merge(indels_patients,indels_controls , by = c("chrom", "snp_position"),all = TRUE)
compare_group_i[is.na(compare_group_i)] <- "no"
compare_group_i$compare[compare_group_i$ALT_control !="no" &  compare_group_i$ALT_patient =="no" ] <- "deletion"
compare_group_i$compare[compare_group_i$ALT_control =="no" &  compare_group_i$ALT_patient !="no" ] <- "addition"
compare_group_i$compare[compare_group_i$ALT_control ==  compare_group_i$ALT_patient ] <- "same"
compare_group_i$compare[compare_group_i$ALT_control !=  compare_group_i$ALT_patient &  compare_group_i$ALT_control !="no" &compare_group_i$ALT_patient !="no"] <- "different"
compare_group_i <- compare_group_i[,c(1,5,11,2,6,12,8,13,14,7,15)]
uncommun_i <- compare_group_i[compare_group_i$compare == "different",]
uncommun_i <- uncommun_i[,-c(2)]
#
deletion_i <- compare_group_i[compare_group_i$compare == "deletion",]
deletion_i <- deletion_i[,-c(2,5)]
deletion_i <- deletion_i[,-c(9)]
deletion_i$ALT_patient <- "."
#
commun <- compare_group[compare_group$compare == "same",]
#filter snps in biomarkers
disease_type_data  <- disease_type_daTa[,c(1,5,7,3)]

names( compare_group_i)[2] <- "gene"
names( compare_group_i)[3] <- "genes"
#disease_type_gene3 <- disease_type_gene[,c(1,5,7)]
annotated_indel5 = sqldf("
  SELECT *
  FROM compare_group_i d1 JOIN disease_type_data d2
  ON d1.gene = d2.gene
  OR d1.genes = d2.gene
")
annotated_indel5 <- annotated_indel5[,-c(2,3)]
annotated_indel5 <- annotated_indel5 %>% 
  mutate_all(funs(str_replace(., "no", ".")))
annotated_indel5<-annotated_indel5 %>% arrange(disease_type)
annotated_indel5 <-annotated_indel5[,c(13,10,11,12,1:9)]


