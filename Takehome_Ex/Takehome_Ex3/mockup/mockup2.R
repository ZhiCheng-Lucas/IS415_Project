# ===================================================================
# Interactive Dual-view Spatial Analysis of Jakarta Educational Infrastructure
# ===================================================================
library(shiny)
library(sf)
library(spdep)
library(ggplot2)
library(dplyr)

# Load pre-processed spatial data
district_metrics <- readRDS("district_metrics.rds")

# Create spatial weights matrix at startup
district_neighbors <- poly2nb(district_metrics, queen = TRUE)
district_weights <- nb2listw(district_neighbors, style = "W")

# ===================================================================
# UI Definition with Dual View
# ===================================================================
ui <- fluidPage(
  titlePanel("Jakarta Schools Spatial Analysis"),
  
  sidebarLayout(
    # Unified sidebar for all controls
    sidebarPanel(
      # Metric selection dropdown
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
      
      # Single significance level control
      sliderInput("significance", "Significance Level",
                  min = 0.01, max = 0.1,
                  value = 0.05, step = 0.01
      ),
      
      # Combined statistical summary
      verbatimTextOutput("analysis_stats")
    ),
    
    # Main panel with two maps side by side
    mainPanel(
      fluidRow(
        column(6,
               h4("LISA (Local Moran's I)", align = "center"),
               plotOutput("lisa_map", height = "500px")
        ),
        column(6,
               h4("Hot Spot (Getis-Ord G*)", align = "center"),
               plotOutput("gstar_map", height = "500px")
        )
      )
    )
  )
)

# ===================================================================
# Server Logic
# ===================================================================
server <- function(input, output, session) {
  
  # Combined reactive calculation for both analyses
  results <- reactive({
    # Extract selected metric
    x <- district_metrics[[input$metric]]
    
    # Create two copies of the data for different analyses
    lisa_results <- district_metrics
    gstar_results <- district_metrics
    
    #------------------------------------------------------------
    # LISA Analysis
    #------------------------------------------------------------
    local_moran <- localmoran(x, district_weights)
    lag_var <- lag.listw(district_weights, x)
    z_var <- scale(x)
    z_lag <- scale(lag_var)
    
    # Initialize LISA classifications
    lisa_clusters <- rep("Not Significant", length(z_var))
    sig_indices <- which(local_moran[, "Pr(z != E(Ii))"] < input$significance)
    
    # Classify LISA clusters
    for(i in sig_indices) {
      if(z_var[i] > 0 && z_lag[i] > 0) lisa_clusters[i] <- "High-High"
      else if(z_var[i] < 0 && z_lag[i] < 0) lisa_clusters[i] <- "Low-Low"
      else if(z_var[i] > 0 && z_lag[i] < 0) lisa_clusters[i] <- "High-Low"
      else if(z_var[i] < 0 && z_lag[i] > 0) lisa_clusters[i] <- "Low-High"
    }
    
    lisa_results$cluster_type <- lisa_clusters
    lisa_results$statistic <- local_moran[, "Ii"]
    lisa_results$p_value <- local_moran[, "Pr(z != E(Ii))"]
    
    #------------------------------------------------------------
    # G* Analysis
    #------------------------------------------------------------
    x_std <- scale(x)
    g_star <- lag.listw(district_weights, x) / scale(x)[,1]
    
    # Calculate z-scores
    n <- length(x)
    W <- listw2mat(district_weights)
    W2 <- W * W
    S1 <- sum(W2)
    mean_g <- sum(W %*% x) / n
    var_g <- S1 * (n - 1) / (n * n - 1)
    z_scores <- (g_star - mean_g) / sqrt(var_g)
    p_values <- 2 * pnorm(-abs(z_scores))
    
    # Classify G* clusters
    gstar_clusters <- rep("Not Significant", length(z_scores))
    gstar_clusters[z_scores > 1.96 & p_values < input$significance] <- "Hot Spot (95%)"
    gstar_clusters[z_scores > 2.58 & p_values < input$significance] <- "Hot Spot (99%)"
    gstar_clusters[z_scores < -1.96 & p_values < input$significance] <- "Cold Spot (95%)"
    gstar_clusters[z_scores < -2.58 & p_values < input$significance] <- "Cold Spot (99%)"
    
    gstar_results$cluster_type <- gstar_clusters
    gstar_results$statistic <- z_scores
    gstar_results$p_value <- p_values
    
    return(list(lisa = lisa_results, gstar = gstar_results))
  })
  
  # LISA Map
  output$lisa_map <- renderPlot({
    res <- results()$lisa
    
    colors <- c(
      "High-High" = "#FF0000",
      "Low-Low" = "#0000FF",
      "High-Low" = "#FF69B4",
      "Low-High" = "#87CEEB",
      "Not Significant" = "#CCCCCC"
    )
    
    ggplot() +
      geom_sf(data = res, aes(fill = cluster_type)) +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      labs(
        fill = "LISA Cluster Type",
        subtitle = paste("Metric:", gsub("_", " ", input$metric))
      )
  })
  
  # G* Map
  output$gstar_map <- renderPlot({
    res <- results()$gstar
    
    colors <- c(
      "Hot Spot (99%)" = "#FF0000",
      "Hot Spot (95%)" = "#FF6666",
      "Cold Spot (99%)" = "#0000FF",
      "Cold Spot (95%)" = "#6666FF",
      "Not Significant" = "#CCCCCC"
    )
    
    ggplot() +
      geom_sf(data = res, aes(fill = cluster_type)) +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      labs(
        fill = "G* Cluster Type",
        subtitle = paste("Metric:", gsub("_", " ", input$metric))
      )
  })
  
  # Combined statistical summary
  output$analysis_stats <- renderPrint({
    res_lisa <- results()$lisa
    res_gstar <- results()$gstar
    
    cat("=== Analysis Summary ===\n\n")
    
    # LISA Statistics
    cat("LISA Cluster Distribution:\n")
    lisa_counts <- table(res_lisa$cluster_type)
    lisa_percentages <- prop.table(lisa_counts) * 100
    for(i in 1:length(lisa_counts)) {
      cat(sprintf("%s: %.1f%% (%d districts)\n",
                  names(lisa_counts)[i],
                  lisa_percentages[i],
                  lisa_counts[i]))
    }
    
    # Global Moran's I
    global_moran <- moran.test(res_lisa[[input$metric]], district_weights)
    cat("\nGlobal Moran's I:\n")
    cat(sprintf("I = %.3f (p-value = %.3f)\n\n",
                global_moran$estimate[1],
                global_moran$p.value))
    
    # G* Statistics
    cat("G* Cluster Distribution:\n")
    gstar_counts <- table(res_gstar$cluster_type)
    gstar_percentages <- prop.table(gstar_counts) * 100
    for(i in 1:length(gstar_counts)) {
      cat(sprintf("%s: %.1f%% (%d districts)\n",
                  names(gstar_counts)[i],
                  gstar_percentages[i],
                  gstar_counts[i]))
    }
    
    cat(sprintf("\nG* Summary:\n"))
    cat(sprintf("Mean z-score = %.3f\n", mean(res_gstar$statistic)))
    cat(sprintf("Z-score range: %.3f to %.3f\n",
                min(res_gstar$statistic), max(res_gstar$statistic)))
  })
}

# Launch the application
shinyApp(ui = ui, server = server)