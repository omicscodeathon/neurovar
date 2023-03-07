
names(annotation1)[9] <- "gene"
expression_gene_interest2_plus_ontology <- merge(expression_gene_interest2, annotation1,
                                                 by.x = "gene", by.y = "gene")
annotated_gene_list <- expression_gene_interest2_plus_ontology[,c(1,2,3,4,5,6,7,11,16,18,19,23,25)]
annotated_gene_list_exp <-unique( annotated_gene_list[,-c(12,13)])

target_gene_annot <-  annotated_gene_list_exp %>% filter( annotated_gene_list_exp$gene =="TUBA4A" )
target_gene_annotion <- unique( target_gene_annot[,-c(5,6,7,9,11)])
# add more info / location
target_gene_annotion2 <- data.frame(t(target_gene_annotion[-1]))
colnames(target_gene_annotion2) <- target_gene_annotion[, 1]


############################"
annotated_gene_list_exp <-unique( annotated_gene_list[,-c(12,13)])

target_gene_annot <-  annotated_gene_list_exp %>% filter( annotated_gene_list_exp$gene =="TUBA4A" )
target_gene_annot_transcript <- target_gene_annot[,c(9,11)]

###################"
#ontology
target_gene_GO <-  annotated_gene_list %>% filter( annotated_gene_list$gene =="TUBA4A" )
target_gene_GO <-target_gene_GO[,c(12,13)]
target_gene_GO_BP <-  target_gene_GO %>% filter(target_gene_GO$`GO domain` =="biological_process" )
colnames(target_gene_GO_BP)[2] <- "BP"
target_gene_GO_BP$number <-  row.names(target_gene_GO_BP)
#target_gene_GO_BP <-target_gene_GO_BP[,c(2)]
target_gene_GO_CC <-  target_gene_GO %>% filter(target_gene_GO$`GO domain` =="cellular_component" )
colnames(target_gene_GO_CC)[2] <- "CC"
target_gene_GO_CC$number <-  row.names(target_gene_GO_CC)
#target_gene_GO_CC <-target_gene_GO_CC[,c(2)]
target_gene_GO_MF <-  target_gene_GO %>% filter(target_gene_GO$`GO domain` =="molecular_function" )
colnames(target_gene_GO_MF)[2] <- "MF"
target_gene_GO_MF$number <-  row.names(target_gene_GO_MF)
#target_gene_GO_MF <-target_gene_GO_MF[,c(2)]
#
#put all data frames into list
df_list <- list(target_gene_GO_BP, target_gene_GO_CC, target_gene_GO_MF)

#merge all data frames in list
df_go <- df_list %>% reduce(full_join, by='number')
df_go <-df_go[,c(2,5,7)]
