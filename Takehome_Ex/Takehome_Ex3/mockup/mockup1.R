# ===================================================================
# Interactive Spatial Analysis of Jakarta Educational Infrastructure
# ===================================================================
# Purpose: 
# - Create interactive visualization of school distribution patterns
# - Allow dynamic switching between LISA and G* analysis
# - Enable real-time adjustment of significance levels
# - Provide statistical summaries of spatial patterns
# ===================================================================

# Required libraries for different functionalities
library(shiny)      # Web application framework
library(sf)         # Spatial data handling
library(spdep)      # Spatial statistics calculations
library(ggplot2)    # Visualization
library(dplyr)      # Data manipulation

# Load pre-processed spatial data
# district_metrics contains all our calculated metrics and geometry
district_metrics <- readRDS("district_metrics.rds")

# Create spatial weights matrix at startup
district_neighbors <- poly2nb(district_metrics, queen = TRUE)
district_weights <- nb2listw(district_neighbors, style = "W")

# ===================================================================
# UI Definition
# ===================================================================
ui <- fluidPage(
  titlePanel("Jakarta Schools Spatial Analysis"),
  
  sidebarLayout(
    # Sidebar for user inputs
    sidebarPanel(
      # Analysis type selection
      # - LISA shows both clusters and outliers
      # - G* focuses on hot/cold spot intensity
      radioButtons("analysis_type", "Analysis Type",
                   choices = c(
                     "LISA (Local Moran's I)" = "lisa",
                     "Hot Spot (Getis-Ord G*)" = "gstar"
                   ),
                   selected = "lisa"
      ),
      
      # Metric selection dropdown
      # - Values match column names in district_metrics
      # - Labels formatted for readability
      selectInput("metric", "Select Metric",
                  choices = c(
                    "School Density" = "school_density",
                    "Elementary Density" = "elementary_density",
                    "Junior High Density" = "junior_high_density",
                    "Senior High Density" = "senior_high_density",
                    "Special Ed Density" = "special_ed_density",
                    "Public Ratio" = "public_ratio"
                  ),
                  selected = "school_density"
      ),
      
      # Significance level slider
      # - Range from 0.01 to 0.1 covers common thresholds
      # - Default 0.05 is standard in statistics
      sliderInput("significance", "Significance Level",
                  min = 0.01, max = 0.1,
                  value = 0.05, step = 0.01
      ),
      
      # Statistical summary output
      verbatimTextOutput("analysis_stats")
    ),
    
    # Main panel for map display
    mainPanel(
      plotOutput("analysis_map", height = "600px")
    )
  )
)

# ===================================================================
# Server Logic
# ===================================================================
server <- function(input, output, session) {
  
  # Reactive calculation of spatial statistics
  # Returns modified district_metrics with new columns for visualization
  results <- reactive({
    # Extract selected metric
    x <- district_metrics[[input$metric]]
    
    if(input$analysis_type == "lisa") {
      #------------------------------------------------------------
      # LISA Analysis Section
      #------------------------------------------------------------
      # Calculate Local Moran's I statistics
      local_moran <- localmoran(x, district_weights)
      
      # Get spatially lagged values for classification
      lag_var <- lag.listw(district_weights, x)
      
      # Standardize variables for comparison
      z_var <- scale(x)
      z_lag <- scale(lag_var)
      
      # Initialize cluster classifications
      clusters <- rep("Not Significant", length(z_var))
      
      # Identify significant locations based on user-selected threshold
      sig_indices <- which(local_moran[, "Pr(z != E(Ii))"] < input$significance)
      
      # Classify significant clusters based on quadrant location
      for(i in sig_indices) {
        if(z_var[i] > 0 && z_lag[i] > 0) clusters[i] <- "High-High"
        else if(z_var[i] < 0 && z_lag[i] < 0) clusters[i] <- "Low-Low"
        else if(z_var[i] > 0 && z_lag[i] < 0) clusters[i] <- "High-Low"
        else if(z_var[i] < 0 && z_lag[i] > 0) clusters[i] <- "Low-High"
      }
      
      # Add results to spatial dataframe
      district_metrics$cluster_type <- clusters
      district_metrics$statistic <- local_moran[, "Ii"]
      district_metrics$p_value <- local_moran[, "Pr(z != E(Ii))"]
      
    } else {
      #------------------------------------------------------------
      # Getis-Ord G* Analysis Section
      #------------------------------------------------------------
      # Standardize input variable
      x_std <- scale(x)
      
      # Calculate G* statistic
      g_star <- lag.listw(district_weights, x) / scale(x)[,1]
      
      # Calculate z-scores using matrix operations for efficiency
      n <- length(x)
      W <- listw2mat(district_weights)
      W2 <- W * W
      S1 <- sum(W2)
      mean_g <- sum(W %*% x) / n
      var_g <- S1 * (n - 1) / (n * n - 1)
      z_scores <- (g_star - mean_g) / sqrt(var_g)
      
      # Calculate two-tailed p-values
      p_values <- 2 * pnorm(-abs(z_scores))
      
      # Classify clusters based on significance and z-score thresholds
      clusters <- rep("Not Significant", length(z_scores))
      clusters[z_scores > 1.96 & p_values < input$significance] <- "Hot Spot (95%)"
      clusters[z_scores > 2.58 & p_values < input$significance] <- "Hot Spot (99%)"
      clusters[z_scores < -1.96 & p_values < input$significance] <- "Cold Spot (95%)"
      clusters[z_scores < -2.58 & p_values < input$significance] <- "Cold Spot (99%)"
      
      # Add results to spatial dataframe
      district_metrics$cluster_type <- clusters
      district_metrics$statistic <- z_scores
      district_metrics$p_value <- p_values
    }
    
    return(district_metrics)
  })
  
  # Generate visualization
  output$analysis_map <- renderPlot({
    res <- results()
    
    # Set appropriate color scheme based on analysis type
    if(input$analysis_type == "lisa") {
      colors <- c(
        "High-High" = "#FF0000",    # Hot spots (red)
        "Low-Low" = "#0000FF",      # Cold spots (blue)
        "High-Low" = "#FF69B4",     # High outliers (pink)
        "Low-High" = "#87CEEB",     # Low outliers (light blue)
        "Not Significant" = "#CCCCCC" # Not significant (gray)
      )
      legend_title <- "LISA Cluster Type"
    } else {
      colors <- c(
        "Hot Spot (99%)" = "#FF0000",  # Strong hot spots
        "Hot Spot (95%)" = "#FF6666",  # Moderate hot spots
        "Cold Spot (99%)" = "#0000FF", # Strong cold spots
        "Cold Spot (95%)" = "#6666FF", # Moderate cold spots
        "Not Significant" = "#CCCCCC"  # Not significant
      )
      legend_title <- "G* Cluster Type"
    }
    
    # Create map using ggplot2
    ggplot() +
      geom_sf(data = res, aes(fill = cluster_type)) +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      labs(
        title = paste(
          if(input$analysis_type == "lisa") "LISA Analysis:" else "Hot Spot Analysis:",
          gsub("_", " ", input$metric)
        ),
        subtitle = paste("Significance Level:", input$significance),
        fill = legend_title
      )
  })
  
  # Generate statistical summary
  output$analysis_stats <- renderPrint({
    res <- results()
    
    cat("Analysis Summary:\n\n")
    
    # Calculate and display cluster distribution
    cluster_counts <- table(res$cluster_type)
    cluster_percentages <- prop.table(cluster_counts) * 100
    
    cat("Cluster Distribution:\n")
    for(i in 1:length(cluster_counts)) {
      cat(sprintf("%s: %.1f%% (%d districts)\n",
                  names(cluster_counts)[i],
                  cluster_percentages[i],
                  cluster_counts[i]))
    }
    
    # Display analysis-specific statistics
    if(input$analysis_type == "lisa") {
      # Global Moran's I for overall spatial autocorrelation
      global_moran <- moran.test(res[[input$metric]], district_weights)
      cat("\nGlobal Moran's I:\n")
      cat(sprintf("I = %.3f (p-value = %.3f)\n",
                  global_moran$estimate[1],
                  global_moran$p.value))
    } else {
      # Summary statistics for G*
      cat("\nG* Statistics:\n")
      cat(sprintf("Mean z-score = %.3f\n", mean(res$statistic)))
      cat(sprintf("Range of z-scores: %.3f to %.3f\n",
                  min(res$statistic), max(res$statistic)))
    }
  })
}

# Launch the application
shinyApp(ui = ui, server = server)