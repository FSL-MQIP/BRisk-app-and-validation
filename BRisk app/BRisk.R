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

# Define function to screen for emetic and anthrax risks 
screen_risks <- function(emetic_genes, anthrax_genes) {
  if (is.na(emetic_genes) || is.na(anthrax_genes)) {
    emetic_risk <- "Missing Data"
    anthrax_risk <- "Missing Data"
  } else {
    if (emetic_genes %in% c("4/4()", "3/4()")) {
      emetic_risk <- "Emetic Disease Risk"
    } else if (emetic_genes %in% c("2/4()", "1/4()")) {
      emetic_risk <- "Insufficient Information to Exclude Emetic Disease Risk"
    } else {
      emetic_risk <- "No Evidence of Emetic Disease Risk"
    }
    
    if (anthrax_genes %in% c("3/3()", "2/3()", "1/3()")) {
      anthrax_risk <- "Anthrax Risk"
    } else {
      anthrax_risk <- "No Evidence of Anthrax Risk"
    }
  }
  
  return(list(emetic_risk = emetic_risk, anthrax_risk = anthrax_risk))
}

# Generate database 
# BTyper data
BTyper3_input = read.csv("Btyper3_Results.csv")
colnames(BTyper3_input)[1] <- "Isolate.Name"
gp_input = read.csv("simulation_input.csv")
database = cbind(BTyper3_input,gp_input[,3:7])
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
      numericInput("n0", "Average initial count (mean) in CFU/mL:", value = 100), 
      ), # Numeric input for "initial count"
      
      div(class = "form-group",
      numericInput("d", "Consumer home storage day:", value = 35),
      ), # Numeric input for "storage day"
      
      div(class = "form-group",
      selectInput("foodmatrix",
                  label="Select a food matrix",
                  choices=c("Milk, pasteurized fluid")),
      ), 
      
      div(class = "form-group",
      fileInput("file", "Input BTyper3 result for a detected B cereus isolate"),  # BTyper3 input for a B cereus isolate
      ),
      
      div(class = "form-group",
      submitButton("Submit", icon("refresh"))
        ) 
      ),
    
    mainPanel(
      plotOutput("hist1"),
      htmlOutput("riskOutput"),
      plotOutput("hist2")
    )
  )
)

# Define server
server <- function(input, output) {
  
  # Input BTyper3 result for a B cereus isolate 
  data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath)
    
    df <- df %>% 
      separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
      separate(Adjusted_panC_Group.predicted_species., into = c("panC_Group","predicted_species"), sep = "\\(") %>%
      mutate(ANI = gsub("\\)", "", ANI),
             panC_Group = gsub("\\)", "", panC_Group),
             predicted_species = gsub("\\)", "", predicted_species))
    
    colnames(df)[8] <- "anthrax_genes"
    colnames(df)[9] <- "emetic_genes"
    
    # Input for risk text 
    emetic_genes <- df$emetic_genes
    anthrax_genes <- df$anthrax_genes
    
    # Filter the database input for rows with the same species as the BTyper3 input
    df$species <- trimws(df$species)
    matching_species_df <- subset(database, species == df$species)
    
    # Simulate HTST milk products along the supply chain 
    ## Set seed
    set.seed(1)
    
    ## Randomly assign isolate names to 1000*length(isolates) units of HTST milk products
    ## All isolates from the same species are equally represented
    isolates <- matching_species_df$Isolate.Name
    sampled_isolates <- character()
    for (isolate in isolates) {
      sampled_isolates <- c(sampled_isolates, rep(isolate, 1000))
    }
    sampled_isolates <- sample(sampled_isolates)
    
    ## Set up dataframe for modeling 1000*length(isolates) units of HTST milk products
    n_sim = 1000*length(isolates)
    ModelData = data.frame(unit_id = rep(seq(1,n_sim)))
    ModelData$isolate <- sampled_isolates
    
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
    ModelData$t_H <- rep(input$d, each = n_sim)
    
    ## Model temperature profiles of 1000*length(isolates) units HTST milk 
    env_cond_time <- matrix(c(rep(0,1000*length(isolates)),
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
    
    # Assign serving size (ml) to 1000*length(isolates) units of HTST milk 
    serving.size<-sample(x = c(rep(x = 244,50),rep(245,25),rep(488,20),rep(732,5)),size = 1000*length(isolates),replace = TRUE)
    ModelData$serving.size = serving.size*0.97
    
    ## Generate simulation input 
    ## Assign growth parameters to 1000*length(isolates) units of HTST milk 
    ModelData$index = match(ModelData$isolate, matching_species_df$Isolate.Name)
    ModelData$Q0 = matching_species_df$Q0[ModelData$index]
    ModelData$Nmax = matching_species_df$Nmax[ModelData$index]
    ModelData$b = matching_species_df$b[ModelData$index]
    ModelData$Tmin = matching_species_df$Tmin[ModelData$index]
    ModelData$Clade = matching_species_df$Clade[ModelData$index]
    
    ## Generate N0 from a Poisson distribution 
    set.seed(42)
    N0 = rpois(n = n_sim, lambda = input$n0)
    ModelData$N0 = N0 
    
    ModelData$Topt = sapply(ModelData$Clade, xopt_func)
    ModelData$mu_opt = (ModelData$b*(ModelData$Topt-ModelData$Tmin))^2
    
    # Run simulation
    for (i in 1:nrow(ModelData)){
      my_primary <- list(mu_opt = ModelData$mu_opt[i], Nmax = ModelData$Nmax[i], N0 = ModelData$N0[i], Q0 = ModelData$Q0[i])
      sec_temperature <- list(model = "reducedRatkowsky", xmin = ModelData$Tmin[i], b = ModelData$b[i], clade = ModelData$Clade[i])
      my_secondary <- list(temperature = sec_temperature)
      growth <- predict_dynamic_growth(times = env_cond_time[i,],
                                       env_conditions = tibble(time = env_cond_time[i,],
                                                               temperature = env_cond_temp[i,]),
                                       my_primary,
                                       my_secondary)
      sim <- growth$simulation
      ModelData$conc[i] = tail(sim$logN, 1)
    }
    
    ModelData$realCFU = 10^ModelData$conc
    ModelData$CFU_per_serve = ModelData$realCFU*ModelData$serving.size
    log_CFU_per_serving <- log10(ModelData$CFU_per_serve)
    
    # Return the required data frame
    return(list(df1 =data.frame(ModelData, log_CFU_per_serving),
                df2 = df,
                df3 = cytotoxicity_input))
    })
  
  # Generate a histogram for the distribution of cfu per serving in all HTST milk units
  output$hist1 <- renderPlot({
    req(data())
    df1 <- data()$df1
    df1$color<-ifelse(test = df1$log_CFU_per_serving>=5,yes = "Above 5 log",no = 
                           ifelse(df1$log_CFU_per_serving>=3,yes = "Between 3 and 5 log",no = "Below 3 log"))
    min_value <- min(df1$log_CFU_per_serving)
    max_value <- max(df1$log_CFU_per_serving)
    breaks <- seq(floor(min_value), ceiling(max_value) + 0.1, by = 0.1)
    finalhist<-ggplot(data = df1,aes(x = log_CFU_per_serving))
    finalhist<-finalhist+
      geom_histogram(data = df1,aes(fill=color),binwidth = 0.1, breaks = breaks)+
      scale_fill_manual("B cereus count per serving", 
                        values = c("Above 5 log"="red3","Between 3 and 5 log"="darkorange1","Below 3 log"="springgreen3"))+
      xlab("CFU per Serving (log scale)") +
      ylab("Number of Servings") +
      ggtitle("Distribution of B cereus count (cfu/serving) in HTST milk products") +
      theme_minimal() + 
      theme(plot.title = element_text(size = 24, face = "bold"),       
            axis.title.x = element_text(size = 22),                    
            axis.title.y = element_text(size = 22),  
            axis.text.x = element_text(size = 22),                    
            axis.text.y = element_text(size = 22),
            legend.text = element_text(size = 22),                     
            legend.title = element_text(size = 22, face = "bold"))    
    return(finalhist)
  })
  
  # Generate risk text
  output$riskOutput <- renderUI({
    req(data())
    df2 <- data()$df2
    emetic_genes <- df2$emetic_genes
    anthrax_genes <- df2$anthrax_genes
    risk_result <- screen_risks(emetic_genes, anthrax_genes)
    
    # Generate risk text
    risk_text <- HTML(paste("<div style='position: absolute; top: 800px; left: 20px;'>",
                            "<h2> This is an isolate from phylogenetic <b>", df2$panC_Group,
                            "</b></h2><h2>With <b><span style='color: darkred;'>", risk_result$emetic_risk, 
                            "</span></b></h2><h2> And <b><span style='color: darkred;'>", risk_result$anthrax_risk, 
                            "</div>", sep = ""))
    
    return(risk_text)
  })
  
  # Generate a histogram for the distribution of normalized cytotoxicity for diarrheal risk 
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
      ggtitle(paste("Histogram of Cytotoxicity Values for All Isolates and Phylogenetic", df2$panC_Group, 
                    "\n", 
                    "(Diarrheal Disease Risk Assessment)")) + 
      theme_minimal() +
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
}    

# Run the application 
shinyApp(ui = ui, server = server)
