#requirement
if (!require("shiny")) install.packages("shiny")
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("purrr")) install.packages("purrr")
if (!require("vcfR")) install.packages("vcfR")
if (!require("bslib")) install.packages("bslib")
if (!require("stringr")) install.packages("stringr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("shinydashboard")) install.packages("shinydashboard")
if (!require("DT")) install.packages("DT")
if (!require("vcfR")) install.packages("vcfR")
if (!require("sqldf")) install.packages("sqldf")
if (!require("fs")) install.packages("fs")
library(shiny)
library(dplyr)
library(readr)
library(tidyverse)
library(purrr)
library(vcfR)
library(bslib)
library(stringr)
library(ggplot2)
library(shinydashboard) 
library(DT)
library(sqldf)
library(data.table)
library(fs)

#source data
annotation <- read_csv("internal_data/annotation.txt") 
full_list <- read_csv("internal_data/full_list.csv")
# app ui
disease_type_data <- full_list %>% filter(disease == "Epilepsy")
disease_type_data <-unique(disease_type_data$disease_type)

types_list <-c( disease_type_data)
gene_data <- full_list %>% filter(disease_type == "adult neuronal ceroid lipofuscinosis")
gene_data <- unique(gene_data$gene)
gene_list <-c(gene_data)




########################
gene_info <- full_list %>% filter(
  full_list$disease=="Epilepsy"
  & 
    full_list$disease_type=="adult neuronal ceroid lipofuscinosis"
  &
    full_list$gene=="CTSF"
)
###
gene_trans <- annotation %>% 
filter(gene == "CTSF") %>% 
  dplyr::select(9, 10, 13:16) %>% 
  unique()


####################   tab 2 ################
##############" tab 2#######################
########"read file
library(readr)

  df <- read.csv("demonstration_data/expression.csv", header = TRUE,sep = ",")
  selected_cols <- df[, c("gene", "pvalue","log2FoldChange")]
  names(selected_cols) <- c("gene", "pvalue", "logfc")
  expression <- cbind(df[, !names(df) %in% c("gene", "pvalue","log2FoldChange")], selected_cols)
  expression$expression_profile <- "Not Significant"
  
  expression$expression_profile[expression$logfc  > 0  & expression$pvalue < 0.1 ] <- "Upregulated genes"
  expression$expression_profile[expression$logfc  < 0  & expression$pvalue < 0.1 ] <- "Downregulated genes"
  expression <- expression[,c("gene","pvalue","logfc","expression_profile" )]
  # filter biomarkers
  biom_list <- full_list %>% filter(full_list$`disease`=="Epilepsy")
  biom_list <- biom_list[,c(1)]
  biom_exp = sqldf("
            SELECT *
            FROM  biom_list d1 JOIN expression d2
            ON d1.gene = d2.gene
          ")
  biom_exp <- biom_exp[,-c(1)]
  ggplot(data=expression, aes(x=logfc, y=-log10(pvalue), col=expression_profile)) +
    geom_point() + theme_minimal()+ ggtitle("Volcano Plot")+
    theme(text = element_text(size = 20)) 
  
  
  ################ tab 3
  ti <- c("A>G","G>A","C>T","T>C")
  # Transveersion (Tv)
  tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")
  
  
  
  ###prepare data
  folder_name <- "control"
  folder_path <- "C:/Users/LENOVO/Desktop/NeuroVar/Neurovar juin2024/NeuroVar2024/demonstration_data/snp"
  folder_path <- file.path(folder_path, folder_name)
  file_paths_control <- fs::dir_ls(folder_path)
  #part I : controls
  file_contenets_control <-file_paths_control %>%
    map(function(path){
      fread(path)%>% as.data.frame()
    })
  
      #add each list item into the environment variable as a tibble.
      # and set each list item to an environment variable
      file_contenets_control%>% map(as_tibble) %>% 
        list2env(envir = .GlobalEnv)
      
      # To join all files together into one data frame
      comined_control <-file_contenets_control %>% map(as_tibble) %>%
        purrr::reduce(full_join)
      
      # delete NAs
      # result : keep commun between all samples
      comined_control <- na.omit(comined_control)
      # and set each list item to an environment variable
      file_contenets_control%>% map(as_tibble) %>% 
        list2env(envir = .GlobalEnv)
      
      # To join all files together into one data frame
      comined_control <-file_contenets_control %>% map(as_tibble) %>%
        purrr::reduce(full_join)
      
      # delete NAs
      # result : keep commun between all samples
      comined_control <- na.omit(comined_control)
      ## identify snps
      
      #
      comined_control$nuSub <- paste0(comined_control$REF,">",comined_control$ALT)
      comined_control$TiTv[comined_control$nuSub %in% ti] <- "Ti"
      comined_control$TiTv[comined_control$nuSub %in% tv] <- "Tv"
      #
      comined_control <-dplyr::rename(comined_control,"snp_controls"="nuSub")
      comined_control <-dplyr::rename(comined_control,"TiTv_controls"="TiTv")
      #comined_control <- comined_control[,c("snp_controls","TiTv_controls")]
      ######"patient
      folder_name2 <- "patient"
      folder_path <- "C:/Users/LENOVO/Desktop/NeuroVar/Neurovar juin2024/NeuroVar2024/demonstration_data/snp"
      
      folder_path2 <- file.path(folder_path, folder_name2)
      file_paths_patient <- fs::dir_ls(folder_path2)
      #
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
        purrr::reduce(full_join)
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
      comined_patient
      ###"""""compare
      compare_group <- merge(comined_patient, comined_control, by = c("#CHROM", "POS"), all = TRUE)
      compare_group[is.na(compare_group)] <- "no"
      compare_group$compare[compare_group$snp_controls !="no" &  compare_group$snp_patient =="no" ] <- "addition"
      compare_group$compare[compare_group$snp_controls =="no" &  compare_group$snp_patient !="no" ] <- "deletion"
      compare_group$compare[compare_group$snp_controls ==  compare_group$snp_patient ] <- "same"
      compare_group$compare[compare_group$snp_controls !=  compare_group$snp_patient &  compare_group$snp_controls !="no" &compare_group$snp_patient !="no"] <- "different"
  #########################################################    
      #final table 
      #compare_group_final <- compare_group[c(1,2,3,4,5,14,15,18,26,27,28)]
      compare_group_final <- compare_group[c("#CHROM", "POS", "REF.x","ALT.x","snp_patient","ID.y","REF.y","ALT.y","compare","TiTv_controls","snp_controls","TiTv_patient")]
      
      names(compare_group_final)[1] <- "Chrom"

      chrom_snppos <- compare_group_final[,c(1,2)]
      names(chrom_snppos)[2] <- "snp_position"
      #
      gene_pos <- unique(annotation[,c(5,6,7,9)] )#>>> filet genes
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
      
      ann_snps3 <- ann_snps2[,c("chrom","start","end","gene","snp_position","ID.y","REF.x", "ALT.x","ALT.y","compare","TiTv_controls","TiTv_patient")]
      names(ann_snps3)[1] <- "Chromosome"
      names(ann_snps3)[2] <- "Start" 
      names(ann_snps3)[3] <- "End"
      names(ann_snps3)[4] <- "Gene"
      names(ann_snps3)[5] <- "SNP Position"
      names(ann_snps3)[6] <-"SNP ID"
      names(ann_snps3)[7] <- "Reference Genome Allele"
      names(ann_snps3)[8] <-  "Control's Allele"
      names(ann_snps3)[9] <-"Patient's Allele"
      names(ann_snps3)[10] <-"Comparison"
      names(ann_snps3)[11] <-"TiTv control"
      names(ann_snps3)[12] <-"TiTv Patient"
      ann_snps3[ann_snps3 == "no"] <- "-"
      ann_snps3[ann_snps3 == "same"] <- "Population specific"
      #####
      biom_list <- full_list %>% filter(full_list$`disease`== "Epilepsy")
      biom_list <- biom_list[,c(1)]
      biom_var = sqldf("
            SELECT *
            FROM  biom_list d1 JOIN ann_snps3 d2
            ON d1.gene = d2.gene
           
          ")
      biom_var <- biom_var[,-c(5)]
      names(biom_var)[1] <- "Gene"
      
####################### indels
      ti <- c("A>G","G>A","C>T","T>C")
      # Transveersion (Tv)
      tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")
      
      
      
      ###prepare data
      folder_name <- "control"
      folder_path <- "C:/Users/LENOVO/Desktop/NeuroVar/Neurovar juin2024/NeuroVar2024/demonstration_data/indel"

      folder_path <- file.path(folder_path, folder_name)
      file_paths_control_indel <- fs::dir_ls(folder_path)
      file_contenets_control_indel <- list()
      file_contenets_control_indel <-file_paths_control_indel %>%
        map(function(path){
          fread(path)
        })
      #
      folder_name2 <- "patient"
      folder_path2 <- "C:/Users/LENOVO/Desktop/NeuroVar/Neurovar juin2024/NeuroVar2024/demonstration_data/indel"
      
      folder_path2 <- file.path(folder_path2, folder_name2)
      file_paths_patients_indel <- fs::dir_ls(folder_path2)
      
      file_contenets_patients_indel <- list()
      file_contenets_patients_indel <-file_paths_patients_indel %>%
        map(function(path){
          fread(path)
        })
      #patient
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
      gene_pos <- unique(annotation[,c(5,6,7,9)] )#>>> filet genes
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
      gene_pos <- unique(annotation[,c(5,6,7,9)] )#>>> filet genes
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
      compare_group_i <- compare_group_i[,c(1,9,10,11,2,13,14,8,15)]

      names(compare_group_i)[1] <- "Chromosome"      
      names(compare_group_i)[2] <- "Start"   
      names(compare_group_i)[3] <- "End"   
          names(compare_group_i)[4] <- "Gene"
          names(compare_group_i)[5] <- "Position"  
          names(compare_group_i)[6] <- "Reference Genome Allele"
          names(compare_group_i)[7] <-  "Control's Allele"
          names(compare_group_i)[8] <-"Patient's Allele"
          names(compare_group_i)[9] <-"Comparison" 
          compare_group_i[compare_group_i == "no"] <- "-"
          compare_group_i[compare_group_i == "same"] <- "Population specific"
      # compare_group_i
      # filter biomarkersbiom_list <- full_list %>% filter(full_list$`disease`== "Epilepsy")
      biom_list <- full_list %>% filter(full_list$`disease`=="Epilepsy")
      biom_list <- biom_list[,c(1)]
      biom_var = sqldf("
            SELECT *
            FROM  biom_list d1 JOIN compare_group_i d2
            ON d1.gene = d2.gene
          ")
      biom_var <- biom_var[,-c(3,4,6,7)]
      biom_var