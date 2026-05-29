# Load packages
library(shiny)
library(tidyverse)
library(tibble)
library(EnvStats)         # to load rtri function 
library(truncnorm)        # to load rtruncnorm function
library(jmuOutlier)       # to load rlaplace function
library(formula.tools)    # to load 'rhs'
library(purrr)            # to load 'map'
library(deSolve)          # to load 'ode'
library(rlang)
library(ggplot2)

# Load utility functions
source("UtilityFunctions_dynamic_growth.R")

# Define function to screen for emetic potential
screen_risks <- function(emetic_genes) {
  if (is.na(emetic_genes)) {
    emetic_potential <- "Missing Data"
  } else {
    if (emetic_genes %in% c("4/4()")) {
      emetic_potential <- "Very High Emetic Disease Potential"
    } else if (emetic_genes %in% c("3/4()")) {
      emetic_potential <- "High Emetic Disease Potential"
    } else if (emetic_genes %in% c("2/4()")) {
      emetic_potential <- "Medium Emetic Disease Potential"
    } else if (emetic_genes %in% c("1/4()")) {
      emetic_potential <- "Low Emetic Disease Potential"
    } else {
      emetic_potential <- "Negligible Emetic Disease Potential"
    }
  }
}

# Generate database 
# BTyper data
BTyper3_input = read.csv("Btyper3_Results.csv")
colnames(BTyper3_input)[1] <- "Isolate.Name"
gp_input = read.csv("simulation_input.csv")
database = cbind(BTyper3_input,gp_input[,3:10])
database <- database %>% 
  separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
  mutate(ANI = gsub("\\)", "", ANI))

# Cytotoxicity data
cytotoxicity_input = read.csv("Cytotoxicity_data.csv")
colnames(cytotoxicity_input)[1] <- "Isolate.Name"

# Define ui
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
        /* Increase font size for the entire interface */
        body {
          font-size: 18px;
        }
        
        /* Increase font size for specific elements */
        h1, h2, h3 {
          font-size: 24px;
        }
        
        p {
          font-size: 16px;
        }
        
        /* Set font size for the 'BRisk' title */
        .brisk-title {
          font-size: 30px;
        }
        
        /* Set spacing for the interface elements */
        .sidebar {
          margin-top: 30px;
          margin-bottom: 30px;
        }
        
        .form-group {
          margin-bottom: 30px; /* Increase this value to increase the space */
        }
      ")
    )
  ),

  titlePanel(tags$h1("BRisk", class = "brisk-title")),
  
  sidebarLayout(
    sidebarPanel(
      div(class = "form-group",
      numericInput("n0", "Initial contamination concentration (CFU/mL):", value = 1), 
      ), # Numeric input for "initial contamination concentration"
      
      div(class = "form-group",
          numericInput("volume", "Container size (mL):", value = 1900), 
      ), # Numeric input for "container size"
      
      div(class = "form-group",
      numericInput("d", "Shelf-life day:", value = 21),
      ), # Numeric input for "shelf-life day"
      
      div(class = "form-group",
      selectInput("foodmatrix",
                  label="Select a food matrix",
                  choices=c("Milk, pasteurized fluid")),
      ), 
      
      div(class = "form-group",
      fileInput("file1", "BRiskTyper result"),  # BRiskTyper result
      ),
      
      div(class = "form-group",
          fileInput("file2", "FastANI result"),  # FastANI result
      ),
      
      div(class = "form-group",
      submitButton("Submit", icon("refresh"))
        ) 
      ),
    
    mainPanel(
      fluidRow(
        column(width = 12,plotOutput("hist1"))
      ),
      
      fluidRow(
        column(width = 12,plotOutput("hist2"))
      ),
      
      fluidRow(
        column(width = 12, uiOutput("summary_text"))
      )
  )
 )
)

# Define server
server <- function(input, output) {
  
  # Input BTyper3 result for a B cereus isolate 
  data <- reactive({
    req(input$file1, input$file2)
    df <- read.csv(input$file1$datapath)
    
    df <- df %>% 
      separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
      separate(Adjusted_panC_Group.predicted_species., into = c("panC_Group","predicted_species"), sep = "\\(") %>%
      mutate(ANI = gsub("\\)", "", ANI),
             panC_Group = gsub("\\)", "", panC_Group),
             predicted_species = gsub("\\)", "", predicted_species))
    
    colnames(df)[9] <- "emetic_genes"
    
    # Input BTyper3 result for ANI
    ANI_file <- read.csv(input$file2$datapath, header = FALSE)
    colnames(ANI_file) <- c("query", "reference", "ANI", "matched_genes", "total_genes")
    ANI_file <- ANI_file[, c("reference", "ANI")]
    ANI_file$reference <- sub("^(PS\\d+).*", "\\1", ANI_file$reference)
    ANI_file$num_id <- as.numeric(sub("PS", "", ANI_file$reference))
    ANI_file <- ANI_file[order(ANI_file$num_id), ]
    ANI_file <- ANI_file[, c("reference", "ANI")]
    
    # Input for risk text 
    emetic_genes <- df$emetic_genes
    
    # Filter the database input for rows with the same species as the BTyper3 input
    df$species <- trimws(df$species)
    database$ANI_new = ANI_file$ANI #Adding new ANI
    matching_species_df <- subset(database, species == df$species)
    # Assigning ANI weight 
    matching_species_df$ANI_wght <- matching_species_df$ANI_new / sum(matching_species_df$ANI_new)
    
    # Simulate HTST milk products along the supply chain
    ## Set seed
    set.seed(1)
    
    ## Assign isolate names to 10000 units of HTST milk products
    ## Isolates from the same species are represented by weight determined by ANI
    n_sim = 10000
    matching_species_df$n_units <- round(n_sim * matching_species_df$ANI_wght)
    sampled_isolates <- rep(matching_species_df$Isolate.Name, matching_species_df$n_units)
    if(length(sampled_isolates) != n_sim){
      sampled_isolates <- sample(sampled_isolates, n_sim)
    }
    
    ModelData <- data.frame(
      unit_id = seq_len(n_sim),
      isolate = sampled_isolates
    )
    # Stage 1: facility storage 
    ## (a)  Sample the temperature distribution
    ModelData$T_F <- rep(runif(n_sim,min=3.5,max=4.5)) #uniform distribution
    ## (b) Sample the storage time (in days) distribution
    ModelData$t_F <- rep(runif(n_sim,min=1,max=2)) #uniform distribution
    
    # Stage 2: transport from facility to retail store
    ## (a)  Sample the temperature distribution
    ModelData$T_T <- rep(rtri(n_sim,min=1.7,max=10.0,mode=4.4)) #triangular distribution
    ## (b) Sample the transportation time (in days) distribution
    ModelData$t_T <- rep(rtri(n_sim,min=1,max=10,mode=5))
    
    # Stage 3: storage/display at retail store
    ## (a)  Sample the temperature distribution
    ModelData$T_S <- rep(rtruncnorm(n_sim,a=-1.4,b=5.4,mean=2.3,sd=1.8)) #truncated normal distribution
    ## (b) Sample the storage time (in days) distribution
    ModelData$t_S <- rep(rtruncnorm(n_sim,a=0.042,b=10.0, mean=1.821,sd=3.3)) #truncated normal distribution
    
    ## Stage 4: transportation from retail store to home
    ## (a)  Sample the temperature distribution
    ModelData$T_T2 <- rep(rtruncnorm(n_sim,a=0,b=10,mean=8.5,sd=1.0)) #truncated normal distribution
    ## (b) Sample the transportation time (in days) distribution 
    ModelData$t_T2 <- rep(rtruncnorm(n_sim,a=0.01,b=0.24, mean=0.04,sd=0.02)) #truncated normal distribution
    
    ## Stage 5: home storage 
    ## (a)  Sample the temperature distribution
    temps <- rep(NA, n_sim)
    for (i in 1:n_sim){
      number <- rlaplace(1,m=4.06,s=2.31)
      while (number > 15 | number < -1) {
        number <- rlaplace(1,m=4.06,s=2.31) #truncated laplace distribution 
      }
      temps[i] <- number
    }
    ModelData$T_H <- temps
    
    ## (b) Define t_H as 35 days for all units
    ModelData$t_H <- rep(35, each = n_sim)
    
    ## Model temperature profiles of 1000 units HTST milk 
    env_cond_time <- matrix(c(rep(0,n_sim),
                              ModelData$t_F, 
                              ModelData$t_F+0.001,
                              ModelData$t_F + ModelData$t_T,
                              ModelData$t_F + ModelData$t_T+0.001,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S+0.001,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2+0.001,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2 + ModelData$t_H), ncol = 10)
    
    env_cond_temp <- matrix(c(ModelData$T_F, 
                              ModelData$T_F,
                              ModelData$T_T,
                              ModelData$T_T,
                              ModelData$T_S,
                              ModelData$T_S,
                              ModelData$T_T2,
                              ModelData$T_T2,
                              ModelData$T_H,
                              ModelData$T_H), ncol = 10)
    
    ## Generate simulation input 
    ## Assign growth parameters to 10000 units of HTST milk 
    ModelData$index = match(ModelData$isolate, matching_species_df$Isolate.Name)
    ModelData$mean_LOG10Q0 = matching_species_df$mean_LOG10Q0[ModelData$index]
    ModelData$sd_LOG10Q0 = matching_species_df$sd_LOG10Q0[ModelData$index]
    ModelData <- ModelData %>%
      mutate(LOGQ0 = rnorm(n(),mean = mean_LOG10Q0,sd = sd_LOG10Q0))
    ModelData$Q0 = 10^ModelData$LOGQ0
    ModelData$mean_Nmax = matching_species_df$mean_Nmax[ModelData$index]
    ModelData$sd_Nmax = matching_species_df$sd_Nmax[ModelData$index]
    ModelData <- ModelData %>%
      mutate(LOGNmax = rnorm(n(),
                             mean = mean_Nmax,
                             sd = sd_Nmax))
    ModelData$Nmax = 10^(ModelData$LOGNmax)
    ModelData$b = matching_species_df$b[ModelData$index]
    ModelData$Tmin = matching_species_df$Tmin[ModelData$index]
    ModelData$Clade = matching_species_df$Clade[ModelData$index]
    
    ## Generate N0 from a Poisson distribution 
    set.seed(42)
    N0 = rpois(n = n_sim, lambda = input$n0*input$volume)
    ModelData$N0 = N0/1900 
    
    ModelData$Topt = sapply(ModelData$Clade, xopt_func)
    ModelData$mu_opt = (ModelData$b*(ModelData$Topt-ModelData$Tmin))^2
    
    # Run simulation
    my_times <- seq(0,35)
    num_iterations <- nrow(ModelData)
    all_simulations <- list()
    for (i in 1:num_iterations) {
      my_primary <- list(mu_opt = ModelData$mu_opt[i], Nmax = ModelData$Nmax[i], N0 = ModelData$N0[i], Q0 = ModelData$Q0[i])
      sec_temperature <- list(model = "reducedRatkowsky", xmin = ModelData$Tmin[i], b = ModelData$b[i], clade = ModelData$Clade[i])
      my_secondary <- list(temperature = sec_temperature)
      growth <- predict_dynamic_growth(times = my_times,
                                       env_conditions = tibble(time = env_cond_time[i,],
                                                               temperature = env_cond_temp[i,]),
                                       my_primary,
                                       my_secondary)
      sim <- growth$simulation
      all_simulations[[i]] <- sim
    }
    
    final_conc <- do.call(rbind, all_simulations)
    dat <- final_conc
    dat_end_of_shelf = subset(dat, time == input$d)
    
    # Return the required data frame
    return(list(df1 = dat_end_of_shelf,
                df2 = df,
                df3 = cytotoxicity_input))
    })
  
  # Generate a histogram for the distribution of cfu per serving in all HTST milk units
  output$hist1 <- renderPlot({
    req(data())
    df1 <- data()$df1
    df1<- df1[df1$logN != -Inf, ]
    df1$color<-ifelse(test = df1$logN>=5,yes = "Above 5 log",no = 
    ifelse(df1$logN>=3,yes = "Between 3 and 5 log",no = "Below 3 log"))
    df1$color <- factor(
      df1$color,
      levels = c("Below 3 log", "Between 3 and 5 log", "Above 5 log")
    )
    breaks <- seq(0, 10, by = 1)
    finalhist<-ggplot(data = df1,aes(x = logN))
    finalhist<-finalhist+
      geom_histogram(data = df1,aes(fill=color),binwidth = 0.1, breaks = breaks)+
      scale_fill_manual(name = expression (italic(B~cereus)~"count per ml"), 
                        values = c("Below 3 log"="springgreen3","Between 3 and 5 log"="darkorange1","Above 5 log"="red3"),
                        breaks = c("Below 3 log", "Between 3 and 5 log", "Above 5 log"))+
      xlab("log CFU per ml") +
      ylab("Number of Units (log scale)") +
      ggtitle(expression("Distribution of " * italic("B. cereus") * " count (log CFU per ml) in HTST milk products")) +
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      theme_classic() + 
      theme(plot.title = element_text(size = 24, face = "bold"),       
            axis.title.x = element_text(size = 22),                    
            axis.title.y = element_text(size = 22),  
            axis.text.x = element_text(size = 22),                    
            axis.text.y = element_text(size = 22),
            legend.text = element_text(size = 22),                     
            legend.title = element_text(size = 22, face = "bold"))    
    return(finalhist)
  })
  
  # Generate a histogram for the distribution of normalized cytotoxicity
  output$hist2 <- renderPlot({
    req(data())
    df2 <- data()$df2
    df3 <- data()$df3
    colnames(df3)[colnames(df3) == "Average_Cell_Viability_F"] <- "Normalized_Cytotoxicity"
    df3$panC_Group <- trimws(df3$panC_Group)
    matching_species_df_ct1 <- subset(df3, panC_Group == df2$panC_Group)
    min_value <- min(df3$Normalized_Cytotoxicity)
    max_value <- max(df3$Normalized_Cytotoxicity)
    breaks <- seq(floor(min_value), ceiling(max_value) + 0.05, by = 0.05)
    ggplot() +
      geom_histogram(data = df3, aes(x = Normalized_Cytotoxicity, fill = "All Isolates"), binwidth = 0.05, alpha = 0.8) +
      geom_histogram(data = matching_species_df_ct1, aes(x = Normalized_Cytotoxicity, fill = "Phylogenetic Group"), binwidth = 0.05) +
      xlab("Cytotoxicity Value") +
      ylab("Number of Isolates") +
      ggtitle(paste("Cytotoxicity Distribution of All Isolates and Phylogenetic", df2$panC_Group, 
                    "\n")) + 
      theme_classic() +
      theme(plot.title = element_text(size = 24, face = "bold"),       
            axis.title.x = element_text(size = 22),                    
            axis.title.y = element_text(size = 22),  
            axis.text.x = element_text(size = 22),                    
            axis.text.y = element_text(size = 22),
            legend.text = element_text(size = 22)) +
      scale_fill_manual(values = c("All Isolates" = "lightblue", "Phylogenetic Group" = "yellow"),
                        labels = c("All Isolates", paste("Phylogenetic", df2$panC_Group, sep = " "))) +
      labs(fill = "")
  })
  
  # Generate text output
  output$summary_text <- renderUI({
    req(data())
    df1<- data()$df1
    df2 <- data()$df2
    df3 <- data()$df3
    
    # exposure 
    logN <- df1$logN
    pct_above5 <- sum(logN > 5, na.rm = TRUE)/10000 * 100
    
    exposure <- if (pct_above5 >= 1) {
      "Very High Exposure"
    } else if (any(logN > 5, na.rm = TRUE)) {
      "High Exposure"
    } else if (any(logN > 3 & logN < 5, na.rm = TRUE)) {
      "Medium Exposure"
    } else if (all(logN < 3, na.rm = TRUE)) {
      "Low Exposure"
    } else {
      "Negligible Exposure"
    }
    
    exposure_color <- case_when(
      exposure == "Very High Exposure" ~ "red",
      exposure == "High Exposure" ~ "darkorange",
      exposure == "Medium Exposure" ~ "yellow",
      exposure == "Low Exposure" ~ "green",
      TRUE ~ "lightgreen"
    )
    
    # emetic
    emetic_potential <- screen_risks(df2$emetic_genes)
    emetic_color <- case_when(
      grepl("Very High", emetic_potential) ~ "red",
      grepl("High", emetic_potential) ~ "darkorange",
      grepl("Medium", emetic_potential) ~ "yellow",
      grepl("Low", emetic_potential) ~ "green",
      TRUE ~ "lightgreen"
    )
    
    # diarrheal
    colnames(df3)[colnames(df3) == "Average_Cell_Viability_F"] <- "Normalized_Cytotoxicity"
    df3$panC_Group <- trimws(df3$panC_Group)
    
    matching <- subset(df3, panC_Group == df2$panC_Group)
    ref_ecdf <- ecdf(df3$Normalized_Cytotoxicity)
    
    matching$percentile_ref <- ref_ecdf(matching$Normalized_Cytotoxicity) * 100
    pc = median(matching$percentile_ref)
    
    diarrheal_potential <- case_when(
      pc >= 80 ~ "Very High Diarrheal Disease Potential",
      pc >= 60 ~ "High Diarrheal Disease Potential",
      pc >= 40 ~ "Medium Diarrheal Disease Potential",
      pc >= 20 ~ "Low Diarrheal Disease Potential",
      TRUE ~ "Very Low Diarrheal Disease Potential"
    )
    
    diarrheal_color <- case_when(
      grepl("Very High", diarrheal_potential) ~ "red",
      grepl("High", diarrheal_potential) ~ "darkorange",
      grepl("Medium", diarrheal_potential) ~ "yellow",
      grepl("Low", diarrheal_potential) ~ "green",
      TRUE ~ "lightgreen"
    )
    
    # final risk score
    if (exposure == "Negligible Exposure" ||
        emetic_potential == "Negligible Emetic Disease Potential") {
      emetic_risk <- "Negligible Emetic Disease Risk"
    } else if (emetic_potential == "Low Emetic Disease Potential" &&
               exposure == "Low Exposure") {
      emetic_risk <- "Low Emetic Disease Risk"
    } else if (emetic_potential == "Low Emetic Disease Potential" &&
               (exposure == "Medium Exposure" || exposure == "High Exposure")) {
      emetic_risk <- "Medium Emetic Disease Risk"
    } else if (emetic_potential == "Medium Emetic Disease Potential" &&
               (exposure == "Low Exposure" || exposure == "Medium Exposure")) {
      emetic_risk <- "Medium Emetic Disease Risk"
    } else if (emetic_potential == "High Emetic Disease Potential" &&
               exposure == "Low Exposure") {
      emetic_risk <- "Medium Emetic Disease Risk"
    } else if (exposure == "Very High Exposure" &&
               emetic_potential == "High Emetic Disease Potential") {
      emetic_risk <- "Very High Emetic Disease Risk"
    } else if (emetic_potential == "Very High Emetic Disease Potential" &&
               (exposure == "Medium Exposure" ||
                exposure == "High Exposure" ||
                exposure == "Very High Exposure")) {
      emetic_risk <- "Very High Emetic Disease Risk"
    } else {
      emetic_risk <- "High Emetic Disease Risk"
    }
    
    if (exposure == "Negligible Exposure") {
      diarrheal_risk <- "Negligible Diarrheal Disease Risk"
    } else if (
      (exposure == "Very High Exposure" &&
       diarrheal_potential == "High Diarrheal Disease Potential") ||
      (diarrheal_potential == "Very High Diarrheal Disease Potential" &&
       exposure %in% c("Medium Exposure", "High Exposure", "Very High Exposure"))
    ) {
      diarrheal_risk <- "Very High Diarrheal Disease Risk"
    } else if (
      (diarrheal_potential == "Very Low Diarrheal Disease Potential" &&
       exposure %in% c("Low Exposure", "Medium Exposure")) ||
      (diarrheal_potential == "Low Diarrheal Disease Potential" &&
       exposure == "Low Exposure")
    ) {
      diarrheal_risk <- "Low Diarrheal Disease Risk"
    } else if (
      (exposure == "Very High Exposure" &&
       diarrheal_potential %in% c("Low Diarrheal Disease Potential",
                                  "Medium Diarrheal Disease Potential")) ||
      (exposure == "High Exposure" &&
       diarrheal_potential %in% c("Medium Diarrheal Disease Potential",
                                  "High Diarrheal Disease Potential")) ||
      (exposure == "Medium Exposure" &&
       diarrheal_potential == "High Diarrheal Disease Potential") ||
      (exposure == "Low Exposure" &&
       diarrheal_potential == "Very High Diarrheal Disease Potential")
    ) {
      diarrheal_risk <- "High Diarrheal Disease Risk"
    } else {
      diarrheal_risk <- "Medium Diarrheal Disease Risk"
    }
    
    risk_color_emetic <- if (emetic_risk == "Very High Emetic Disease Risk") {
      "red"
    } else if (emetic_risk == "High Emetic Disease Risk") {
      "darkorange"
    } else if (emetic_risk == "Medium Emetic Disease Risk") {
      "yellow"
    } else if (emetic_risk == "Low Emetic Disease Risk") {
      "green"
    } else {
      "lightgreen"
    }
    
    risk_color_diarrheal <- if (diarrheal_risk == "Very High Diarrheal Disease Risk") {
      "red"
    } else if (diarrheal_risk == "High Diarrheal Disease Risk") {
      "darkorange"
    } else if (diarrheal_risk == "Medium Diarrheal Disease Risk") {
      "yellow"
    } else if (diarrheal_risk == "Low Diarrheal Disease Risk") {
      "green"
    } else {
      "lightgreen"
    }
    
    # print text
    HTML(paste0(
      "<div>",
      
      "<p style='color:black; font-weight:bold; font-size:22px;'>",
      "Summary Report",
      "</p>",
      
      "<p style='color:", exposure_color, "; font-weight:bold; font-size:20px;'>",
      exposure,
      "</p>",
      
      "<p style='color:", emetic_color, "; font-weight:bold; font-size:20px;'>",
      emetic_potential,
      "</p>",
      
      "<p style='color:", diarrheal_color, "; font-weight:bold; font-size:20px;'>",
      diarrheal_potential,
      "</p>",
      
      "<hr>",
      
      "<p style='color:black; font-weight:bold; font-size:20px;'>",
      "Overall Risk Score:",
      "</p>",
      
      "<p style='color:", risk_color_emetic, "; font-weight:bold; font-size:18px;'>",
      emetic_risk,
      "</p>",
      
      "<p style='color:", risk_color_diarrheal, "; font-weight:bold; font-size:18px;'>",
      diarrheal_risk,
      "</p>",
      
      "</div>"
    ))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
