### indel visualization
### requirements
library(bedr)
library(tidyverse)
library(ggplot2)

####
indel_data <- read.vcf("demo_data/indel/patients/p2_SRR12878718.vcf") 
indels <- indel_data$vcf
alt_length <- nchar(indels$ALT)
ref_length <- nchar(indels$REF)
indel_length <- alt_length - ref_length
indel_signs <- sign(indel_length)
count_table <- data.frame(length = as.factor(indel_length), variant = ifelse(indel_signs == -1, "Deletion", "Insertion")) %>% 
  table() %>% 
  as.data.frame() %>% 
  setNames(c("length", "variant", "count"))
## plot
ggplot(count_table, aes(x = length, y = count, fill = variant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("pink", "cyan"), labels = c("Deletion", "Insertion")) +
  xlab("Indel length") +
  ylab("Count")

#################
one_indels <- indels[,c(1,2,4,5)]
one_indels <-one_indels %>%  mutate_all(funs(str_replace(., "chr", "")))

# add controls
one_indels2 = sqldf("
  SELECT *
  FROM one_indels d1 JOIN annotated_snps2 d2
  ON d1.CHROM = d2.chrom
  AND d1.POS = d2.POS_snp_control
  AND d1.REF = d2.REF_control
")
one_indels2 <- one_indels2[,c(1, 8, 6, 7, 2, 10, 3, 4, 12)]
#####
#type

one_indels2$length_REF<-nchar(gsub("[^A-Z]","",one_indels2$REF))
one_indels2$length_ALT<-nchar(gsub("[^A-Z]","",one_indels2$ALT))
one_indels2$length_ALTc<-nchar(gsub("[^A-Z]","",one_indels2$ALT_control))
# the ref = input
# vs control !!
one_indels2$type_vs_control[one_indels2$length_ALT < one_indels2$length_ALTc ] <- "deletion"
one_indels2$type_vs_control[one_indels2$length_ALT > one_indels2$length_ALTc ] <- "addition"
one_indels2$type_vs_control[one_indels2$length_ALT == one_indels2$length_ALTc 
                            & one_indels2$ALT != one_indels2$ALT_control ] <- "polymorphism"
one_indels2$type_vs_control[one_indels2$ALT == one_indels2$ALT_control ] <- "identical"

#vs ref
one_indels2$type_vs_reference[one_indels2$length_ALT < one_indels2$length_REF ] <- "deletion"
one_indels2$type_vs_reference[one_indels2$length_ALT > one_indels2$length_REF ] <- "addition"
one_indels2$type_vs_reference[one_indels2$length_ALT == one_indels2$length_REF
                            & one_indels2$ALT != one_indels2$REF ] <- "polymorphism"
one_indels2$type_vs_reference[one_indels2$ALT == one_indels2$REF ] <- "identical"

# remove identical # depend input
#one_indels2_d <- one_indels2 %>%
 # filter( one_indels2$type_vs_control !="identical")
# filter biomarker
one_indels3= sqldf("
  SELECT *
  FROM one_indels2 d1 JOIN disease_type_data d2
  ON d1.gene = d2.gene
")
one_indels3 <- one_indels3[,-c(10, 11, 12, 15, 16, 17)]

one_indels3<-one_indels3 %>% arrange(disease_type)


