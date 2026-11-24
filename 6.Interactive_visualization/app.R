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
df_merged <- inner_join(df_counts_t, df_activity, by = "sample_id")

# MEA variable names
#n_activity <- ncol(df_activity_numeric) - 1
#mea_cols <- (ncol(df_merged) - n_activity + 1):ncol(df_merged)
#mea_names <- colnames(df_merged)[mea_cols]
# Identify MEA variable names correctly (numeric columns only)
mea_names <- df_activity %>%
  select(where(is.numeric)) %>%
  colnames()
# ----------------------------
# Shiny UI
# ----------------------------
ui <- fluidPage(
  titlePanel("Normalized gene expression vs activity"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("gene", "Enter gene name:", value = ""),
      selectInput("mea_var", "Select MEA variable:", choices = mea_names)
    ),
    
    mainPanel(
      plotOutput("scatterPlot"),
      verbatimTextOutput("corrMsg"),
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
      mea = df_merged[[mea]],
     cell_line  = factor(df_merged$cell_line, levels=c("ADNP_1","YY1_1",
							"EHMT1_4_KS","ARID1B_3",
							"SMARCB1_1","SMARCB1_2",
							"SMARCB1_1_KSS","SMARCB1_2_KSS","CHD8_1")),   
     color = df_merged$color   # color per sample
    )
  })
  
  output$scatterPlot <- renderPlot({
    data <- selected_data()
    if(is.null(data)) return(NULL)
 
  # Create a mapping between Gene and color (unique pairs)
    color_map <- unique(data[, c("cell_line", "color")])
    cell_colors <- setNames(color_map$color, color_map$cell_line)
   
    ggplot(data, aes(x = mea, y = gene)) +
      geom_point(aes(fill = cell_line),shape = 21,size = 5, color = "black") +             # bigger black dots
     geom_smooth(method = "lm", col = "black", linetype = "dashed", se = FALSE, inherit.aes = FALSE,
                aes(x = mea, y = gene)) +  # regression line across all data
      scale_fill_manual(values = cell_colors, name = "Cell line") +  
      xlab(paste0("MEA variable: ", input$mea_var)) +
      ylab(paste0("Gene: ", input$gene)) +
      theme_minimal(base_size = 14) +
      theme(
       legend.title = element_text(size = 13, face = "bold"),
       legend.text = element_text(size = 11),
       legend.position = "right"
    )
  })
  
  output$corrMsg <- renderPrint({
    data <- selected_data()
    if(is.null(data)) return(NULL)
    
    cor_test <- cor.test(data$gene, data$mea, method = "spearman")
    cat(sprintf("Spearman correlation: %.3f, p-value: %.3g", cor_test$estimate, cor_test$p.value))
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
