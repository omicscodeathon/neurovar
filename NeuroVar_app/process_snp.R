###prepare data

#part I : controls
#path
#file_paths_control <- fs::dir_ls("data/snp/control")
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



