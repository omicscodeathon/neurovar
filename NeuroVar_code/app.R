
#source("requirements.R")
#source("disease_gene_code.R")
#source("expression.R")
#source("discover_gene.R")
#source("indel.R")
#source("snp.R")
#source("one_sample_indel.R")
##################
ui <- bootstrapPage(
  navbarPage(theme = shinytheme("simplex"),
             collapsible = TRUE,
             HTML('<a style="text-decoration:none;
               cursor:default;
                    color:#FC2E20;
                    " class="active" href="#">NeuroVar</a>'),
             id="nav",
             windowTitle ="NeuroVar",
             tabPanel(h2("Biomarker"),
                      sidebarLayout(
                        sidebarPanel(
                          pickerInput(inputId ="disease_n",
                                      label = "Select the disease of interest:",
                                      choices = c("disease1", "N/A"),
                                      # choices = c(unique(full_list$disease), "N/A"),
                                      selected = "N/A"
                          ),
                          br(),
                          
                          span(shiny::tags$i(h6("disease Type")),
                               pickerInput(inputId ="disease_t",
                                           label = "Select the disease Type:",
                                           choices = c(unique(types_list$disease_type), "N/A"),
                                           selected = "N/A"
                               )
                          ),
                          selectInput(inputId = "target_gene2",
                                      label = "Gene of interest:",
                                      choices =c(disease_type_gene$gene, "N/A"),
                                      selected ="N/A"
                          )
                          
                        ),
                        mainPanel(
                          fluidRow(
                            h1("The following genes are associated with the disease"),
                            plotOutput(outputId = "karyotype")  
                          ),
                          br(),br(),
                          fluidRow(
                            box("about",h2("About the gene"),width=6,
                                DT::dataTableOutput("gene_info")
                                # h2("The gene expression profile :"),
                                #textOutput("exp_profile") 
                                
                            ),
                            box("trans",width=6, h2("gene's transcript"),
                                DT::dataTableOutput("gene_trans")
                            )#,
                            #box("onto",width=6, h2("The gene ontologies :"),
                            #    DT::dataTableOutput("gene_GO")
                            #)
                          )
                          
                        ))
                      
             ),#tab1
             tabPanel(h2("Control vs patient"),
                      fluidRow(h1("Expression"),
                               box(width = 4, 
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
                                   )
                               ),
                               box(
                                 width = 8,
                                 plotOutput(outputId = "volcano_plot")
                                 
                               )
                               
                      ),# exp FR
                      fluidRow(h1("SNPs"),
                               box("table1",DTOutput("table_snps"))
                      ), # snp FR
                      fluidRow(h1("Indels"),
                               box("table2", DT::dataTableOutput("table_indels")) 
                      ) # indel FR
                      
             ),# tab2
             tabPanel(h2("One Sample"),
                      
                      fluidRow(
                        box(" ")
                      ), # snp FR
                      fluidRow("Indels in Biomarkers", 
                               box(DT::dataTableOutput("table_one_indels")) 
                      ) # indel FR
                      
                      
             )# tab3
  ))






server <- function(input, output, session) {
  

  ###############################################################################################################
  # all genes 
  output$karyotype <- renderPlot({
    gene.symbols <- c(disease_type_gene$gene)
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                             filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
    seqlevelsStyle(genes) <- "UCSC"
    
    #viz
    
   # kp <- plotKaryotype(genome="hg38")
   # kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "horizontal",
   #               r1=0.5, cex=0.8, adjust.label.position = FALSE)
  })
  
  ###################################tab 1###################################################################
  output$exp_table <- DT::renderDataTable({
    exp_tabl <-  expression_gene_interest[,c(1,18,15,21)]
  })
  #volcano plot
  output$volcano_plot <- renderPlot({
    ggplot(data=expression_gene_interest, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
  })
  
  ############################## tab2 ########################################################################

 #table
  output$gene_info <- DT::renderDataTable({

    target_gene_annotion2

  })
  output$gene_trans <- DT::renderDataTable({

    target_gene_annot_transcript
  })

  #
  output$gene_GO <- DT::renderDataTable({
   df_go
  })

  output$table_all_genes_onto <- DT::renderDataTable({
    annotated_gene_list_GO <-annotated_gene_list[,c(1,12,13)]
  })


  
  output$exp_profile <- renderText({
    if (gene_expression$pvalue < 0.05){
      expression_profile <- "The gene is differentially expressed"
    } else{
      expression_profile <- "The gene is NOT differentially expressed"
    }
  })
  ######################################### tab3 snps ########################################################

  output$table_snps <- renderDT(annotated_snps4,
                           filter = "top",
                           class = "display nowrap compact" # style
  )
 ######################################### tab 4 indels###################################################################
  #source("process_indels.R")
 output$table_indels <- renderDT({
   annotated_snps5
   })
  output$table_one_indels <- renderDT({
    one_indels3
  })

}

# Run the application
shinyApp(ui = ui, server = server)




