# Project Title: "Analyzing Gene Expression Time Course Data in Yeast: A Shiny App Project"

# Introduction:
# In this project, we will explore a time course gene expression dataset from Paul Spellman's lab at Stanford. 
# This microarray dataset contains transcript levels from the yeast Saccharomyces cerevisiae genome, 
# which respond to various cyclins, reflecting different stages of the cell cycle. Specifically, we will focus 
# on the cdc15 time course experiment and conduct various data analysis and visualization tasks using R. 
# We will also create a Shiny web application to interactively explore gene expression patterns across different time points.

# Data Loading and Exploration:
# We start by loading the dataset into R, which contains 6,178 genes and 77 time points or arrays.
# We isolate the cdc15 experiment samples for further analysis.

# Load the required libraries
library(ggplot2)
library(shiny)

# Load the data
data <- read.table("C:/Users/sweth/OneDrive/Documents/R learning/Lab2/spellman/spellman.txt",
                   header = TRUE, row.names = 1)

# Isolate only the cdc15 experiment (samples 23-46)
cdc15_data <- data[, 23:46]

# Correlation Matrix Visualization:
# We calculate the correlation matrix between time points using Pearson's correlation.
# A correlation plot is generated using ggplot2, with labeled axes and a color legend.

# Calculate the correlation matrix while avoiding NA values
dat.cor <- cor(cdc15_data, use = 'pairwise.complete.obs')

# Create a layout for the plot
layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2), 5, 2, byrow = TRUE)
       
       # Set outer margins to a smaller value
       par(oma = c(3, 0, 3, 0))
       
       # Define a custom color palette
       cx <- rev(colorRampPalette(c("red", "white", "blue"))(25))
       
       # Create a sequence of values for the legend
       leg <- seq(min(dat.cor, na.rm = TRUE), max(dat.cor, na.rm = TRUE), length.out = 10)
       
       # Create the correlation plot
       image(dat.cor, main = "Correlation plot Expression/Stages of life cycle",
             axes = FALSE, col = cx)
       
       # Add x and y axes with labels
       axis(1, at = seq(0, 1, length.out = ncol(dat.cor)), labels = dimnames(dat.cor)[[2]],
            cex.axis = 0.9, las = 2)
       
       axis(2, at = seq(0, 1, length.out = ncol(dat.cor)), labels = dimnames(dat.cor)[[2]],
            cex.axis = 0.9, las = 2)
       
       # Create the color legend with a smaller plot region
       # Adjust the plot region as needed
       par(plt = c(0.1, 0.5, 0.1, 0.5))
       image(as.matrix(leg), col = cx, axes = FALSE, xlab = 'Pearson correlation')
       tmp <- round(leg, 2)
       axis(1, at = seq(0, 1, length.out = length(leg)), labels = tmp, cex.axis = 1)
       
       # Imputing Missing Values:
       # We select a specific gene, YAL002W (VPS8), and impute missing values with the row mean. 
       # This gene is involved in endosomal vesicle tethering and fusion.
       
       # Select the gene YAL002W (VPS8) and assign the missing values with the row mean
       gene_to_impute1 <- as.numeric(cdc15_data["YAL002W", ])
       gene_to_impute1[is.na(gene_to_impute1)] <- mean(gene_to_impute1, na.rm = TRUE)
       
       # Check the mean
       gene_mean <- mean(gene_to_impute1, na.rm = TRUE)
       
       # Gene Expression Profile Plot:
       # We create a profile plot of gene YAL002W (VPS8) over time. The x-axis shows time points without the "cdc15_" prefix.
       
       # Extract column names from cdc15_data
       sample_names <- colnames(cdc15_data)
       time_points <- as.numeric(sub("cdc15_", "", sample_names))
       x_labels <- seq(10, max(time_points), by = 20)
       
       # Create a data frame with time_points and gene_to_impute
       profile_data <- data.frame(Time = time_points, Expression = gene_to_impute1)
       
       # Create the profile plot using ggplot2
       ggplot(data = profile_data, aes(x = Time, y = Expression)) +
         geom_line(linewidth = 1) +
         labs(
           title = "Expression Profile of YAL002W (VPS8) Over Time",
           x = "Time Points",
           y = "Expression Level"
         ) +
         scale_x_continuous(breaks = x_labels)
       
       # Shiny Web Application:
       # We create a Shiny web application to allow users to select and correlate any time points across all genes.
       # Users can choose X and Y variables, as well as point colors for the scatter plot.
       
       # Define the UI
       ui <- fluidPage(
         sidebarLayout(
           sidebarPanel(
             selectInput('xcol', 'X Variable', colnames(cdc15_data)),
             selectInput('ycol', 'Y Variable', colnames(cdc15_data)),
             selectInput('point_color', 'Point Color',
                         choices = c("Red" = "red", "Blue" = "blue", "Green" = "green",
                                     "Yellow" ="yellow","Pink"="pink","Orange"="orange",
                                     "Purple"="purple","Black" = "black"))
           ),
           mainPanel(
             plotOutput('scatterplot')
           )
         )
       )
       
       # Define the Server
       server <- function(input, output) {
         selectedData <- reactive({
           cdc15_data[, c(input$xcol, input$ycol, input$colorcol)]
         })
         output$scatterplot <- renderPlot({
           data <- selectedData()
           if (!is.null(data)) {
             plot(data[, 1], data[, 2],
                  xlab = input$xcol,
                  ylab = input$ycol,
                  col = input$point_color,
                  pch = 19,
                  main = paste("Scatter Plot:", input$xcol, "vs", input$ycol))
           }
         })
       }
       
       # Running the Shiny app
       shinyApp(ui = ui, server = server)
       