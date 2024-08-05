
UI <- shinyUI({
  navbarPage(
    theme = bs_theme(version = 3, bg = "white", fg = "black", primary ="#ff4f00"),
    collapsible = TRUE,
    HTML('<a style="text-decoration:none;
               cursor:default;
                    color:#ff4f00;
                    " class="active" href="#">NeuroVar</a>'),
    id="nav",
    windowTitle ="NeuroVar",
    
    tabPanel("Biomarker",
             sidebarLayout(
               sidebarPanel(
                 selectInput(inputId ="disease_n",
                             label = "Select the disease of interest:",
                             choices = c(unique(full_list$disease), "N/A"),
                             selected = "N/A"
                 ),
                 
                 selectInput(inputId ="disease_t",
                             label = "Select the disease Type:",
                             choices = "N/A",
                             selected = "N/A"
                 ),
                 selectInput(inputId = "target_gene",
                             label = "Gene of interest:",
                             choices = "N/A",
                             selected ="N/A"
                 )
               ),
               mainPanel(
                 box(h2("About the gene"),width=12,
                     DT::dataTableOutput("gene_infos")
                 ),
                 box(width=12, h2("gene's transcript"),
                     DT::dataTableOutput("gene_transs")
                     
                 )
               )
             )
    ),
    tabPanel("Expression",
             sidebarLayout(
               sidebarPanel(
                 # import file
                 fileInput('target_upload', 'Choose file to upload',
                           accept = c(
                             'text/csv',
                             'text/comma-separated-values',
                             '.csv'
                           )),
                 radioButtons("separator","Separator: ",choices = c(";",",",":"), selected=";",inline=TRUE),
                 #actionButton("submit", label = "Submit"),        
                 # define gene column
                 # Select input controls for column selection
                 selectInput("col1", "Select Gene Column", ""),
                 selectInput("col2", "Select P-value Column", ""),
                 selectInput("col3", "Select LogFC Column", ""),
                 
                 # filters p val log fc
                 h3("Note: Make Sure the disease is selected in the firs tab ! "),
                 
                 # filters p val log fc
                 h4("Define the adjusted p-value and LogFC value to identify the differentially expressed genes"),
                 sliderInput("pval",
                             "Adjusted P_value:",
                             min = 0,
                             max = 0.05,
                             value=0.01
                 ),
                 br(),
                 br(),
                 sliderInput("log",
                             "Log value:",
                             min = 0,
                             max = 5,
                             value=2
                 )
               ),
               mainPanel(
                 #expression table
                 DT::dataTableOutput("expression"),
                 # volcano
                 plotOutput(outputId = "volcano_plot")
                 
                 
               )
             )),
    tabPanel("Variants",
             sidebarLayout(
               sidebarPanel(
                 textInput("folder_path", label = "Enter folder path:"),
                 "Note: Make sure the path contain two folders named 'control' and 'patient'",
                 radioButtons("vtype","variant type: ",choices = c("SNP","Indels"), selected="SNP"),
                 actionButton("submit_button", "Submit"),
               ),
               mainPanel(
                 #"comparison",
                 DTOutput("comined")
               )
             )
    )
    
  )

})

