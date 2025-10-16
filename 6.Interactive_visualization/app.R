###############################################
# app.R
# Gene vs MEA Variable Scatter Plot
# Standalone version (run in R or deploy to shinyapps.io)
###############################################

# ----------------------------
# Load libraries
# ----------------------------
library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# ----------------------------
# Data URLs
# ----------------------------
counts_url <- "https://raw.githubusercontent.com/sofiapuvogelvittini/MEA_seq_extended/master/0.Data/RNAseq_data/control_normalized_neuronal_counts/PRPS_voom_norm_control.txt"
activity_url <- "https://raw.githubusercontent.com/sofiapuvogelvittini/MEA_seq_extended/master/0.Data/MEA_data/Processed/mutants_normalized/patients_normalized_hier_PC13.csv"

# ----------------------------
# Load counts
# ----------------------------
df_counts <- read.table(
  counts_url,
  header = TRUE,
  sep = "",  # whitespace-separated
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# rownames are genes
genes <- rownames(df_counts)

# clean sample names (column names)
colnames(df_counts) <- gsub("\\.", "-", colnames(df_counts))
colnames(df_counts) <- str_replace_all(colnames(df_counts), "SMARCB1_KSS_CRISPR", "SMARCB1_KSS-CRISPR")

# transpose: rows = samples, columns = genes
df_counts_t <- as.data.frame(t(df_counts))
colnames(df_counts_t) <- genes
df_counts_t$sample_id <- rownames(df_counts_t)
rownames(df_counts_t) <- NULL

# ----------------------------
# Load activity
# ----------------------------
df_activity <- read_csv(activity_url, show_col_types = FALSE)
df_activity <- as.data.frame(df_activity)

# remove last column if cluster
df_activity <- df_activity[-ncol(df_activity)]

# only numeric MEA variables + sample_id
df_activity$sample_id <- df_activity$sample
df_activity_numeric <- df_activity %>% select_if(is.numeric)
df_activity_numeric$sample_id <- df_activity$sample_id

# ----------------------------
# Merge counts and activity
# ----------------------------
df_merged <- inner_join(df_counts_t, df_activity_numeric, by = "sample_id")

# MEA variable names
n_activity <- ncol(df_activity_numeric) - 1
mea_cols <- (ncol(df_merged) - n_activity + 1):ncol(df_merged)
mea_names <- colnames(df_merged)[mea_cols]

# ----------------------------
# Shiny UI
# ----------------------------
ui <- fluidPage(
  titlePanel("Gene vs MEA Variable Scatter Plot"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("gene", "Enter gene name:", value = ""),
      selectInput("mea_var", "Select MEA variable:", choices = mea_names)
    ),
    
    mainPanel(
      plotOutput("scatterPlot"),
      verbatimTextOutput("warningMsg")
    )
  )
)

# ----------------------------
# Shiny server
# ----------------------------
server <- function(input, output, session) {
  
  selected_data <- reactive({
    req(input$gene, input$mea_var)
    gene <- input$gene
    mea <- input$mea_var
    
    if(!gene %in% colnames(df_merged)) return(NULL)
    
    data.frame(
      gene = df_merged[[gene]],
      mea = df_merged[[mea]]
    )
  })
  
  output$scatterPlot <- renderPlot({
    data <- selected_data()
    if(is.null(data)) return(NULL)
    
    ggplot(data, aes(x = mea, y = gene)) +
      geom_point(color = "blue") +
      geom_smooth(method = "lm", col = "red", se = FALSE) +
      xlab(paste0("MEA variable: ", input$mea_var)) +
      ylab(paste0("Gene: ", input$gene)) +
      theme_minimal()
  })
  
  output$warningMsg <- renderPrint({
    if(!input$gene %in% colnames(df_merged)){
      cat("Warning: gene not found in dataset!")
    }
  })
}

# ----------------------------
# Run the app
# ----------------------------
shinyApp(ui, server)
