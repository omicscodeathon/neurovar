### indel visualization
### requirements
library(bedr)
library(tidyverse)
library(ggplot2)

####
snp_data <- read.vcf("demo_data/snp/patient/p1_SRR12878707.vcf") 
snp_data <- snp_data$vcf
#################
snp_data  <-snp_data [,c(1:5)]
snp_data  <- snp_data  %>%  mutate_all(funs(str_replace(., "chr", "")))
# ti tv

# add controls

#####
#type  snp_data
snp_data$type <- 'unknown'
compare_group[is.na(compare_group)] <- "no"
compare_group$compare[compare_group$snp_controls !="no" &  compare_group$snp_patient =="no" ] <- "deletion"
compare_group$compare[compare_group$snp_controls =="no" &  compare_group$snp_patient !="no" ] <- "addition"
compare_group$compare[compare_group$snp_controls ==  compare_group$snp_patient ] <- "identical"
compare_group$compare[compare_group$snp_controls !=  compare_group$snp_patient 
                      &  compare_group$snp_controls !="no" 
                      &compare_group$snp_patient !="no"] <- "polymorphism"

# the ref = input
# vs control !!

#vs ref

# remove identical # depend input
one_indels2_d <- one_indels2 %>%
  filter( one_indels2$type_vs_control !="identical")
# filter biomarker
one_indels3= sqldf("
  SELECT *
  FROM one_indels2 d1 JOIN disease_type_data d2
  ON d1.gene = d2.gene
")
one_indels3 <- one_indels3[,-c(10, 11, 12, 15, 16, 17)]

one_indels3<-one_indels3 %>% arrange(disease_type)


