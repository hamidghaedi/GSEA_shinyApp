# GSEA_shinyApp
An application to do GSEA from  RNA-seq raw count

```R
library(shiny)
library(DESeq2)
library(tibble)
library(dplyr)

# Load data
rna <- read.csv("Uromol1_CountData.v1.csv", header = T, sep = ",", row.names = 1)
uromol_clin <- read.csv("uromol_clinic.csv", sep = ",", header = T, row.names = 1)

# Ensure that row names/column names in both dataframes are the same 
uromol_clin <- uromol_clin[colnames(rna),]


# Set grouping variable
uromol_clin$isTa <- as.factor(ifelse(uromol_clin$Stage == "Ta", "Ta", "non_Ta"))
levels(uromol_clin$isTa)

# Create UI elements
ui <- fluidPage(
  titlePanel("DE Analysis App"),
  sidebarLayout(
    sidebarPanel(
      fileInput("rna", "Select RNA file"),
      fileInput("uromol_clin", "Select clinical data file")
    ),
    mainPanel(
      verbatimTextOutput("summary")
    )
  )
)

# Create server logic
server <- function(input, output) {
  # Set maximum file size to 300MB
  options(shiny.maxRequestSize = 300 * 1024^2)
  
  # Reactive expression for data input
  data_input <- reactive({
    rna <- read.csv(input$rna$datapath, header = T, sep = ",", row.names = 1)
    uromol_clin <- read.csv(input$uromol_clin$datapath, sep = ",", header = T, row.names = 1)
    # Ensure that row names/column names in both dataframes are the same 
    uromol_clin <- uromol_clin[colnames(rna),]
    
    # Set grouping variable
    uromol_clin$isTa <- as.factor(ifelse(uromol_clin$Stage == "Ta", "Ta", "non_Ta"))
    
    # Return data frames
    list(rna = rna, uromol_clin = uromol_clin)
  })
  
  # Reactive expression for DE analysis
  de_analysis <- reactive({
    rna <- data_input()$rna
    uromol_clin <- data_input()$uromol_clin
    
    # Make DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(countData = rna,
                                  colData = uromol_clin,
                                  design = ~ isTa)
    # Perform DE analysis
    dds <- DESeq(dds)
    
    # Extract results
    res <- results(dds)
    res <- cbind(as.data.frame(rowData(dds)), res)
    res <- res[, -duplicated(colnames(res))]
    res <- res %>% 
      filter(baseMean > 1) %>%
      mutate(log2FoldChange = log2FoldChange,
             adj_pval = padj)
    
    # Return DE results
    res
  })
  
  # Output DE analysis results
  output$summary <- renderPrint({
    de_analysis()
  })
}

# Run the app
shinyApp(ui, server)


```
