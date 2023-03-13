###prepare data

#part I : controls
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
#uncommun <- compare_group[compare_group$compare == "different",]
#addition <- compare_group[compare_group$compare == "addition",]
#deletion <- compare_group[compare_group$compare == "deletion",]
#commun <- compare_group[compare_group$compare == "same",]

# dataframe to shown 
#uncommun_simple <- data.frame(uncommun$`#CHROM`, uncommun$POS, uncommun$ID.x, 
#                              uncommun$REF.x, uncommun$ALT.x, uncommun$ALT.y)
#rename columns
#names(uncommun_simple) <- c('CHROM','POS','ID','REF','ALT_patient', 'ALT_control')

#final table 
compare_group_final <- compare_group[c(1,2,3,4,5,14,15,18,26,27,28)]
names(compare_group_final)[1] <- "Chrom"


##############################################################################################""
chrom_snppos <- compare_group_final[,c(1,2)]
names(chrom_snppos)[2] <- "snp_position"
names(chrom_snppos)[1] <- "chrom"
#
gene_pos <- unique(annotation1[,c(5,6,7,9)] )#>>> filet genes
names(gene_pos)[1] <- "chrom"
names(gene_pos)[2] <- "start"
names(gene_pos)[3] <- "end"
#•
chrom_snppos  <- chrom_snppos  %>%
  mutate_all(funs(str_replace(., "chr", "")))
compare_group_final  <-compare_group_final %>%
  mutate_all(funs(str_replace(., "chr", "")))
#
ann_snps = sqldf("
  SELECT *
  FROM gene_pos d1 JOIN chrom_snppos d2
  ON d1.chrom = d2.chrom
  AND d2.snp_position < d1.end
  AND d2.snp_position > d1.start
")
ann_snps <- ann_snps[,-c(5)]
#merge with snp info
ann_snps2 = sqldf("
  SELECT *
  FROM ann_snps d1 JOIN compare_group_final d2
  ON d1.chrom = d2.chrom
  AND d1.snp_position = d2.POS
")

ann_snps3 <- ann_snps2[,c(1,4,5,8,9,13,10,16)]
names(ann_snps3)[4] <- "SNP ID"
names(ann_snps3)[5] <- "Reference Genome Allele"
names(ann_snps3)[6] <-  "Control's Allele"
names(ann_snps3)[7] <-"Patient's Allele"
names(ann_snps3)[2] <- "gene"
#filter snps in biomarkers
disease_type_data  <- disease_type_daTa[,c(1,5,7,3)]


#disease_type_gene3 <- disease_type_gene[,c(1,5,7)]
ann_snps4 = sqldf("
  SELECT *
  FROM ann_snps3 d1 JOIN disease_type_data d2
  ON d1.gene = d2.gene
")
ann_snps4<- ann_snps4[,-c(9)]
ann_snps4<-  ann_snps4%>% 
  mutate_all(funs(str_replace(., "no", ".")))
ann_snps4<- ann_snps4 %>% arrange(disease_type)
ann_snps4 <- ann_snps4[,c(11,2,1,9,10,3:8)]

##########################
#regions
names(annotation2)[15] <- "gene"

#######
ann_snps__region = sqldf("
  SELECT *
  FROM ann_snps4 d1 JOIN annotation2 d2
  ON d1.gene = d2.gene
")

#
ann_snps__region[,c(27)]<-"region"
names(ann_snps__region)[27] <- "region"
ann_snps__region[,c(27)]<-"."
# to remove
names(ann_snps__region)[12]<-"coding_start"
names(ann_snps__region)[13]<-"coding_end"
names(ann_snps__region)[20]<-"exon_end"
names(ann_snps__region)[21]<-"exon_start"
names(ann_snps__region)[22]<-"5utr_start"
names(ann_snps__region)[23]<-"5utr_end"
names(ann_snps__region)[24]<-"3utr_start"
names(ann_snps__region)[25]<-"3utr_end"
#
ann_snps__region$region[ann_snps__region$snp_position  <ann_snps__region$coding_end
                      & ann_snps__region$snp_position  > ann_snps__region$coding_start] <- "coding_region"
ann_snps__region$region[ann_snps__region$snp_position  <ann_snps__region$exon_end
                        & ann_snps__region$snp_position  > ann_snps__region$exon_start] <- "exon_region"
ann_snps__region$region[ann_snps__region$snp_position  <ann_snps__region$`5utr_end`
                        & ann_snps__region$snp_position  > ann_snps__region$`5utr_start`] <- "exon_region"
ann_snps__region$region[ann_snps__region$snp_position  <ann_snps__region$`3utr_end`
                        & ann_snps__region$snp_position  > ann_snps__region$`3utr_start`] <- "exon_region"
ann_snps__region$region[ann_snps__region$snp_position  <ann_snps__region$`CDS end`
                        & ann_snps__region$snp_position  > ann_snps__region$`CDS start`] <- "CDS_region"
ann_snps__region$region[ann_snps__region$snp_position  <ann_snps__region$`cDNA coding end`
                        & ann_snps__region$snp_position  > ann_snps__region$`cDNA coding start`] <- "CDS_region"

#filter 
final_region  <- ann_snps__region[,-c(12:26)]
final_region <- unique(final_region)
# rmove . but only if duplicated

