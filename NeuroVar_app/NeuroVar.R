#source("requirements.R")
##################
ui <- bootstrapPage(
  navbarPage(theme = shinytheme("flatly"),
             collapsible = TRUE,
             HTML('<a style="text-decoration:none;
               cursor:default;
                    color:#FFFFFF;
                    " class="active" href="#">app name</a>'),
             id="nav",
             windowTitle ="app name",
             
             sidebarLayout(
               sidebarPanel(
                 pickerInput(inputId ="disease_n",
                             label = "Select the disease of interest:",
                             choices = c(unique(full_list$disease), "N/A"),
                             selected = "N/A"
                 ),
                 br(),
                 #
                 span(shiny::tags$i(h6("disease Type")),
                 pickerInput(inputId ="disease_t",
                             label = "Select the disease Type:",
                             choices = c(unique(types_list$disease_type), "N/A"),
                             selected = "N/A"
                 )
                 )
               ),
               mainPanel(
                 h3("The following genes are associated with the disease"),
                 plotOutput(outputId = "karyotype"),
                 
                 h1("Expression profile of the genes of interest"),
                 
                 span(shiny::tags$i(h6("Define the p-value to identify the differentially expressed genes")), style="color:#045a8d"),
                 br(),
                 sliderInput("pval",
                             "P_value:",
                             min = 0,
                             max = 1,
                             value=0.05
                 ),
                 br(),
                 span(shiny::tags$i(h6("Define the LogFC to identify the differentially expressed genes")), style="color:#045a8d"),
                 br(),
                 sliderInput("log",
                             "Log value:",
                             min = 0,
                             max = 5,
                             value=2
                 ),
                 plotOutput(outputId = "volcano_plot"),
                 br(),
                 br(),
                 br(),
                 DT::dataTableOutput("table_all_genes_exp"),
                 br(),
                 br(),
                 br(),
                 h1("Ontologies of the genes of interest"),
                 DT::dataTableOutput("table_all_genes_onto"),
                 
                 ############
                 #
                 h1("Filter One gene data"),
                 h3("The following genes are associated with the disease"),
                 
                 selectInput(inputId = "target_gene",
                             label = "Gene of interest:",
                             choices =c(disease_type_gene$gene, "N/A"),
                             selected ="N/A"
                 ),
                 # output  expression profile
                 h2("The gene expression profile :"),
                 textOutput("exp_profile"),
                 br(),
                 br(),
                 br(),
                 # output snps
                 h2("SNPs in target gene :"),
                 
                 pickerInput(inputId ="snp_region",
                                  label = "Filter SNPs in  a Genomic region:",
                                  choices = c("Genomic Coding region","CDS","cDNA coding", "Exon","5'UTR","3'UTR"),
                                  selected = "N/A"),
                 tableOutput("snp_gene_table"),
                 br(),
                 br(),
                 br(),
                 # output indels
                 h2("Indelss in target gene :"),
                 br(),
                 br(),
                 br(),
                 # output cnvs
                 h2("CNVs in target gene :"),
                 br(),
                 br(),
                 br(),
               )
             
             
          )
))




server <- function(input, output, session) {
  # filter disease + type
  #disease
  #full_list <- read_excel("gene-disease/full_list.xlsx")
  full_list<- full_list[-1,]
  names(full_list)[3] <- "disease_type"
  names(full_list)[1] <- "gene"
  names(full_list)[10] <- "disease"
  
  #disease subtype

  disease_type_daTa <- full_list %>% filter(full_list$disease == "Epilepsy")
  types_list<- as.data.frame( unique(disease_type_daTa$disease_type))
  names(types_list)[1] <- "disease_type"
  # gene
  disease_type_gene <-  disease_type_daTa %>% filter(disease_type_daTa$disease_type == "epilepsy")
  ###############################################################################################################
  # all genes 
  output$karyotype <- renderPlot({
    gene.symbols <- c(disease_type_gene$gene)
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                             filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
    seqlevelsStyle(genes) <- "UCSC"
    
    #viz
    
    kp <- plotKaryotype(genome="hg38")
    kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "horizontal",
                  r1=0.5, cex=0.8, adjust.label.position = FALSE)
    })
  ###############################################################################################################
  # expression
  #expression <- read_csv("data/expression.csv")
  names(expression)[3] <- "gene"
  expression$diffexpressed <- "Not Significant"

  #expression$diffexpressed[expression$log2FoldChange  > input$log  & expression$pvalue < input$pval ] <- "Upregulated genes"
  #expression$diffexpressed[expression$log2FoldChange  < -(input$log)  & expression$pvalue < input$pval ] <- "Downregulated genes"
  expression$diffexpressed[expression$log2FoldChange  > 2  & expression$pvalue < 0.05 ] <- "Upregulated genes"
  expression$diffexpressed[expression$log2FoldChange  < -2  & expression$pvalue < 0.05 ] <- "Downregulated genes"
  #
  #filter target genes
  expression_gene_interest <-merge(disease_type_gene, expression, by.x = "gene", by.y = "gene")
  expression_gene_interest2 <- expression_gene_interest[,c(1,5,7,8,15,18,21)]
  ##
  names(annotation1)[9] <- "gene"
  expression_gene_interest2_plus_ontology <- merge(expression_gene_interest2, annotation1,
                                                   by.x = "gene", by.y = "gene")
  annotated_gene_list <- expression_gene_interest2_plus_ontology[,c(1,2,3,4,5,6,7,11,16,18,19,23,25)]
  
  ###############################################################################################################
  #table
  output$table_all_genes_exp <- DT::renderDataTable({
    annotated_gene_list_exp <-unique( annotated_gene_list[,-c(12,13)])
  })
  output$table_all_genes_onto <- DT::renderDataTable({
    annotated_gene_list_GO <-annotated_gene_list[,c(1,12,13)]
  })
  #volcano plot
  output$volcano_plot <- renderPlot({
  ggplot(data=expression_gene_interest, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
  })
  ##############################################################################################################
  #source("process_snp.R")
    ##############################################################################################################
  ######################################### one gene ########################################################
  ##############################################################################################################

  input_gene = "CASR" #input
  #filter annotate gene
  
  #annotation1 <- read_csv("annotation/annotation1.txt")
  gene_annotation <-  annotation1 %>% filter(annotation1$gene == input_gene) 
  genes_location <-unique( gene_annotation[c(5,6,7,9)])
  #extract start and end 
  start_gene <- genes_location$`Gene start (bp)`
  end_gene <- genes_location$`Gene end (bp)`
  chromosome_gene <- genes_location$`Chromosome/scaffold name`

  # determine expression profile
  gene_expression <- expression %>% filter(expression$gene == input_gene)
  output$exp_profile <- renderText({
    if (gene_expression$pvalue < 0.05){
      expression_profile <- "The gene is differentially expressed"
    } else{
      expression_profile <- "The gene is NOT differentially expressed"
    }
  })
  ##############################################################################################################
  
  # filter SNPs in that gene
  # filter snp in gene
  snp_gene_filter <-compare_group_final %>% filter(compare_group_final$Chrom==chromosome_gene 
                                                   & compare_group_final$POS < end_gene
                                                   & compare_group_final$POS > start_gene)
  

  #output
  output$snp_gene_table <- renderTable({
    
    # real output 
    #snp_gene_filter
    
    #example table 
    example_snp_table <- compare_group_final[1:10, ]
    
    # To add add count rows
    # if count == 0 >> render text no snps in the gene
  })
  ##############################################################################################################
  
  
  }

# Run the application
shinyApp(ui = ui, server = server)




