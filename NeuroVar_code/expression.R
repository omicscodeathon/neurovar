
# expression


expression$diffexpressed <- "Not Significant"

#expression$diffexpressed[expression$log2FoldChange  > input$log  & expression$pvalue < input$pval ] <- "Upregulated genes"
#expression$diffexpressed[expression$log2FoldChange  < -(input$log)  & expression$pvalue < input$pval ] <- "Downregulated genes"
expression$diffexpressed[expression$log2FoldChange  > 2  & expression$pvalue < 0.05 ] <- "Upregulated genes"
expression$diffexpressed[expression$log2FoldChange  < -2  & expression$pvalue < 0.05 ] <- "Downregulated genes"
#
#filter target genes
expression_gene_interest <-merge(disease_type_gene, expression, by.x = "gene", by.y = "gene")
expression_gene_interest2 <- expression_gene_interest[,c(1,5,7,8,15,18,21)]


#############"""

#filter gene

gene_annotation <-  annotation1 %>% filter(annotation1$`Gene name` == input_gene) 
genes_location <-unique( gene_annotation[c(5,6,7,9)])
#extract start and end 
start_gene <- genes_location$`Gene start (bp)`
end_gene <- genes_location$`Gene end (bp)`
chromosome_gene <- genes_location$`Chromosome/scaffold name`

# determine expression profile
gene_expression <- expression %>% filter(expression$gene == input_gene)
