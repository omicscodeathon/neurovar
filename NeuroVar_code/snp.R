###prepare data

#part I : controls
#path

#loop
file_contenets_control <- list()
file_contenets_control <-file_paths_control %>%
  map(function(path){
    fread(path)
  })
#add each list item into the environment variable as a tibble.
# and set each list item to an environment variable
file_contenets_control%>% map(as_tibble) %>% 
  list2env(envir = .GlobalEnv)
# To join all files together into one data frame7
comined_control <-file_contenets_control %>% map(as_tibble) %>%
  reduce(full_join)
# delete NAs
# result : keep commun between all samples
comined_control <- na.omit(comined_control)
## identify snps
# Transition (Ti)
ti <- c("A>G","G>A","C>T","T>C")
# Transveersion (Tv)
tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")
#
comined_control$nuSub <- paste0(comined_control$REF,">",comined_control$ALT)
comined_control$TiTv[comined_control$nuSub %in% ti] <- "Ti"
comined_control$TiTv[comined_control$nuSub %in% tv] <- "Tv"
#
comined_control <-dplyr::rename(comined_control,"snp_controls"="nuSub")
comined_control <-dplyr::rename(comined_control,"TiTv_controls"="TiTv")
### Part II : patients
#path
#file_paths_patient <- fs::dir_ls("data/snp/patient")
file_contenets_patient <- list()
file_contenets_patient <-file_paths_patient %>%
  map(function(path){
    fread(path)
  })
#add each list item into the environment variable as a tibble.
# and set each list item to an environment variable
file_contenets_patient%>% map(as_tibble) %>% 
  list2env(envir = .GlobalEnv)
# To join all files together into one data frame7
comined_patient <-file_contenets_patient %>% map(as_tibble) %>%
  reduce(full_join)
# delete NAs
# result : keep commun between all samples
comined_patient <- na.omit(comined_patient)
## identify snps
comined_patient$nuSub <- paste0(comined_patient$REF,">",comined_patient$ALT)
comined_patient$TiTv[comined_patient$nuSub %in% ti] <- "Ti"
comined_patient$TiTv[comined_patient$nuSub %in% tv] <- "Tv"
#
comined_patient<-dplyr::rename(comined_patient,"snp_patient"="nuSub")
comined_patient<-dplyr::rename(comined_patient,"TiTv_patient"="TiTv")

# Part III: compare groups
compare_group <- merge(comined_patient, comined_control, by = c("#CHROM", "POS"),all = TRUE)
compare_group[is.na(compare_group)] <- "no"
compare_group$compare[compare_group$snp_controls !="no" &  compare_group$snp_patient =="no" ] <- "deletion"
compare_group$compare[compare_group$snp_controls =="no" &  compare_group$snp_patient !="no" ] <- "addition"
compare_group$compare[compare_group$snp_controls ==  compare_group$snp_patient ] <- "same"
compare_group$compare[compare_group$snp_controls !=  compare_group$snp_patient &  compare_group$snp_controls !="no" &compare_group$snp_patient !="no"] <- "different"
uncommun <- compare_group[compare_group$compare == "different",]
addition <- compare_group[compare_group$compare == "addition",]
deletion <- compare_group[compare_group$compare == "deletion",]
commun <- compare_group[compare_group$compare == "same",]

# dataframe to shown 
uncommun_simple <- data.frame(uncommun$`#CHROM`, uncommun$POS, uncommun$ID.x, 
                              uncommun$REF.x, uncommun$ALT.x, uncommun$ALT.y)
#rename columns
names(uncommun_simple) <- c('CHROM','POS','ID','REF','ALT_patient', 'ALT_control')
#final table 
compare_group_final <- compare_group[c(1,2,3,4,5,14,15,18,26,27,28)]
names(compare_group_final)[1] <- "Chrom"







############"
##############################################################################################""
chrom_snppos <- compare_group_final[,c(1,2)]
names(chrom_snppos)[2] <- "snp_position"
names(chrom_snppos)[1] <- "chrom"
#
gene_pos <- unique(annotation1[,c(5,6,7,9)] )#>>> filet genes
names(gene_pos)[1] <- "chrom"
names(gene_pos)[2] <- "start"
names(gene_pos)[3] <- "end"


annotated_snps = sqldf("
  SELECT *
  FROM gene_pos d1 JOIN chrom_snppos d2
  ON d1.chrom = d2.chrom
  AND d2.snp_position < d1.end
  AND d2.snp_position > d1.start
")
annotated_snps <- annotated_snps[,-c(5)]
#merge with snp info
annotated_snps2 = sqldf("
  SELECT *
  FROM annotated_snps d1 JOIN compare_group_final d2
  ON d1.chrom = d2.chrom
  AND d1.snp_position = d2.POS
")

annotated_snps3 <- annotated_snps2[,c(1,4,5,8,9,13,10,16)]
names(annotated_snps3)[4] <- "SNP ID"
names(annotated_snps3)[5] <- "Reference Genome Allele"
names(annotated_snps3)[6] <-  "Control's Allele"
names(annotated_snps3)[7] <-"Patient's Allele"

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


