# ===================================================================
# Jakarta Educational Infrastructure Spatial Analysis
# ===================================================================
# Purpose: 
#   This application enables interactive spatial analysis of educational 
#   infrastructure distribution in Jakarta. It provides two types of analyses:
#   1. LISA (Local Indicators of Spatial Association) for identifying clusters
#      and spatial outliers
#   2. Getis-Ord G* for identifying hot and cold spots
#
# Key Features:
#   - Interactive visualization of school distribution patterns
#   - Dynamic switching between LISA and G* analysis methods
#   - Real-time adjustment of significance levels
#   - Comprehensive statistical summaries
#   - Multiple educational metrics analysis
#
# Data Requirements:
#   - Input data (district_metrics.rds) must contain:
#     * Spatial geometry for each district
#     * Density metrics for different school types
#     * Public-private ratio data
#
# Output:
#   - Interactive map showing spatial clusters/patterns
#   - Statistical summaries of spatial relationships
#   - Cluster distribution analysis
# ===================================================================

# Required libraries
library(shiny)      # Web application framework for interactive features
library(sf)         # Spatial data handling and manipulation
library(spdep)      # Spatial statistics calculations
library(ggplot2)    # Advanced data visualization
library(dplyr)      # Data manipulation and transformation

# Load and prepare spatial data
# district_metrics.rds should be a spatial dataframe containing:
# - Geometry column for district boundaries
# - Metrics columns for various school densities
# - Additional educational infrastructure metrics
district_metrics <- readRDS("district_metrics.rds")

# Create spatial weights matrix
# - Uses queen contiguity (shared borders or vertices)
# - Converts to listw object for spatial calculations
district_neighbors <- poly2nb(district_metrics, queen = TRUE)
district_weights <- nb2listw(district_neighbors, style = "W")

# ===================================================================
# UI Definition
# ===================================================================
ui <- fluidPage(
  titlePanel("Jakarta Schools Spatial Analysis"),
  
  sidebarLayout(
    # Sidebar Panel: Contains all user input controls
    sidebarPanel(
      # Analysis Type Selection
      # Options:
      # - LISA: Identifies clusters and spatial outliers
      # - G*: Identifies hot and cold spots
      radioButtons("analysis_type", "Analysis Type",
                   choices = c(
                     "LISA (Local Moran's I)" = "lisa",
                     "Hot Spot (Getis-Ord G*)" = "gstar"
                   ),
                   selected = "lisa"
      ),
      
      # Metric Selection
      # Different density and ratio metrics for analysis
      selectInput("metric", "Select Metric",
                  choices = c(
                    "School Density" = "school_density",         # All schools
                    "Elementary Density" = "elementary_density", # Elementary schools
                    "Junior High Density" = "junior_high_density", # Junior high
                    "Senior High Density" = "senior_high_density", # Senior high
                    "Special Ed Density" = "special_ed_density",   # Special education
                    "Public Ratio" = "public_ratio"               # Public vs private
                  ),
                  selected = "school_density"
      ),
      
      # Significance Level Control
      # - Range: 0.01 to 0.1
      # - Default: 0.05 (standard statistical significance)
      sliderInput("significance", "Significance Level",
                  min = 0.01, max = 0.1,
                  value = 0.05, step = 0.01
      ),
      
      helpText("Adjust parameters above to update the analysis.")
    ),
    
    # Main Panel: Contains visualization and statistics
    mainPanel(
      # Map Output
      # Height set to 600px for good visibility of spatial patterns
      plotOutput("analysis_map", height = "600px"),
      
      # Summary Statistics Panel
      wellPanel(
        fluidRow(
          # Cluster Distribution Statistics
          column(6,
                 h5("Cluster Distribution", align = "center"),
                 verbatimTextOutput("cluster_stats")
          ),
          # Analysis-Specific Metrics
          column(6,
                 h5("Analysis Metrics", align = "center"),
                 verbatimTextOutput("specific_stats")
          )
        )
      )
    )
  )
)

# ===================================================================
# Server Logic
# ===================================================================
server <- function(input, output, session) {
  
  # Reactive Spatial Analysis Calculator
  # Purpose: Performs spatial analysis based on user inputs
  # Returns: Modified district_metrics with added analysis columns
  results <- reactive({
    # Extract selected metric for analysis
    x <- district_metrics[[input$metric]]
    
    if(input$analysis_type == "lisa") {
      #------------------------------------------------------------
      # LISA (Local Moran's I) Analysis
      #------------------------------------------------------------
      # Calculate Local Moran's I statistics
      local_moran <- localmoran(x, district_weights)
      
      # Calculate spatially lagged values
      lag_var <- lag.listw(district_weights, x)
      
      # Standardize variables for comparison
      z_var <- scale(x)          # Standardize original values
      z_lag <- scale(lag_var)    # Standardize lagged values
      
      # Initialize cluster classifications
      clusters <- rep("Not Significant", length(z_var))
      
      # Identify significant locations
      sig_indices <- which(local_moran[, "Pr(z != E(Ii))"] < input$significance)
      
      # Classify clusters based on quadrant location
      # - High-High: High value surrounded by high values
      # - Low-Low: Low value surrounded by low values
      # - High-Low: High value surrounded by low values
      # - Low-High: Low value surrounded by high values
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
      # Getis-Ord G* Analysis
      #------------------------------------------------------------
      # Standardize input variable
      x_std <- scale(x)
      
      # Calculate G* statistic
      g_star <- lag.listw(district_weights, x) / scale(x)[,1]
      
      # Calculate z-scores using matrix operations
      n <- length(x)
      W <- listw2mat(district_weights)    # Convert weights to matrix
      W2 <- W * W                         # Square weights matrix
      S1 <- sum(W2)                       # Sum of squared weights
      mean_g <- sum(W %*% x) / n          # Mean G* statistic
      var_g <- S1 * (n - 1) / (n * n - 1) # Variance of G*
      z_scores <- (g_star - mean_g) / sqrt(var_g)
      
      # Calculate two-tailed p-values
      p_values <- 2 * pnorm(-abs(z_scores))
      
      # Classify clusters based on significance and z-score thresholds
      # Uses standard normal distribution critical values
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
  
  # Map Visualization Output
  # Creates the spatial visualization using ggplot2
  output$analysis_map <- renderPlot({
    res <- results()
    
    # Define color schemes for different analysis types
    colors <- if(input$analysis_type == "lisa") {
      c(
        "High-High" = "#FF0000",    # Hot spots (red)
        "Low-Low" = "#0000FF",      # Cold spots (blue)
        "High-Low" = "#FF69B4",     # High outliers (pink)
        "Low-High" = "#87CEEB",     # Low outliers (light blue)
        "Not Significant" = "#CCCCCC" # Not significant (gray)
      )
    } else {
      c(
        "Hot Spot (99%)" = "#FF0000",  # Strong hot spots
        "Hot Spot (95%)" = "#FF6666",  # Moderate hot spots
        "Cold Spot (99%)" = "#0000FF", # Strong cold spots
        "Cold Spot (95%)" = "#6666FF", # Moderate cold spots
        "Not Significant" = "#CCCCCC"  # Not significant
      )
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
        fill = if(input$analysis_type == "lisa") "LISA Cluster Type" else "G* Cluster Type"
      )
  })
  
  # Cluster Distribution Output
  # Displays the distribution of clusters and their percentages
  output$cluster_stats <- renderPrint({
    res <- results()
    
    cat("Cluster Distribution:\n")
    cat("----------------\n")
    cluster_counts <- table(res$cluster_type)
    cluster_percentages <- prop.table(cluster_counts) * 100
    
    for(i in 1:length(cluster_counts)) {
      cat(sprintf("%s:\n  %.1f%% (%d districts)\n",
                  names(cluster_counts)[i],
                  cluster_percentages[i],
                  cluster_counts[i]))
    }
  })
  
  # Analysis-Specific Statistics Output
  # Displays different statistics based on analysis type
  output$specific_stats <- renderPrint({
    res <- results()
    
    if(input$analysis_type == "lisa") {
      # Calculate and display Global Moran's I
      global_moran <- moran.test(res[[input$metric]], district_weights)
      cat("LISA Statistics:\n")
      cat("---------------\n")
      cat(sprintf("Global Moran's I: %.3f\n", global_moran$estimate[1]))
      cat(sprintf("P-value: %.3f\n", global_moran$p.value))
    } else {
      # Display G* summary statistics
      cat("G* Statistics:\n")
      cat("--------------\n")
      cat(sprintf("Mean z-score: %.3f\n", mean(res$statistic)))
      cat(sprintf("Z-score range:\n  Min: %.3f\n  Max: %.3f\n",
                  min(res$statistic), max(res$statistic)))
    }
  })
}

# Launch the application
shinyApp(ui = ui, server = server)