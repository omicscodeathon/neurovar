
#disease


full_list<- full_list[-1,]
names(full_list)[3] <- "disease_type"
names(full_list)[1] <- "gene"
names(full_list)[10] <- "disease"

#disease subtype
#library(dplyr)
disease_type_daTa <- full_list %>% filter(full_list$disease == "Epilepsy")
types_list<- as.data.frame( unique(disease_type_daTa$disease_type))
names(types_list)[1] <- "disease_type"
# gene
disease_type_gene <-  disease_type_daTa %>% filter(disease_type_daTa$disease_type == "epilepsy")



