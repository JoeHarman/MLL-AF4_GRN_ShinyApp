#ui.R

##### Load packages ----------------------------------------------------------------------------------------

library(shiny)
library(shinythemes)
library(DT)
library(knitr)
library(plotly)
library(rfigshare)
library(tidyverse)
library(igraph)
library(ggnetwork)
library(viridis)
library(intergraph)



##### User Interface  ----------------------------------------------------------------------------------------

shinyUI(
  
  navbarPage("MLL-AF4 SEM GRN", fluid = T,
    theme = shinytheme("cosmo"),
    
    
    ##### TAB - description --------------------------------------------------------------------------------------------
    tabPanel("Description",
      
      titlePanel("Interactive navigation of MLL-AF4 SEM GRN paper"),
     
      ##### Main panel ---------------------------------------------------------------------------------------------
      mainPanel(
        h3("The purpose of this application"),
        "This application is supplementary to Harman et al, 2020, providing an interactive user interface to explore transcriptomic data and gene regulatory networks of MLL-AF4 leukemia in SEM cells.",
        br(),
        includeMarkdown(knitr::knit("GRN_paper_affiliations.Rmd"))
      )
    ),
    ##### TAB -  RNA-seq --------------------------------------------------------------------------------------------
    tabPanel(
      "RNA-seq",
      
      titlePanel("RNA-seq expression changes in response to drug and siRNA perturbations"),
      
      ##### Side bar -------------------------------------------------------------------------------------------
      sidebarLayout(
        sidebarPanel(
          position = "left", width=3,
          selectizeInput('Gene_Name_RNA', 
                         label = "Select gene:",
                         choices = NULL)
        ),
        
        ##### Main panel ---------------------------------------------------------------------------------------------
        mainPanel(
          plotOutput("RNA_bp", width = "90%"),
          br(),
          DT::dataTableOutput('table_RNA')
        )
      )
    ),
    
    ##### TAB -  Network --------------------------------------------------------------------------------------------
    tabPanel(
      "GRN navigation",
      
      titlePanel("Predicted direct upstream regulators and downstream targets"),
      
      ##### Side bar -------------------------------------------------------------------------------------------
      sidebarLayout(
        sidebarPanel(
          position = "left",width=3,
          selectizeInput('Gene_Name_GRN', 
                         label = "Select node:",
                         choices = NULL),
          radioButtons(
            "up_down", "Show up or downstream targets: ", 
            choices = c(
              "Downstream targets" = "down",
              "Upstream regulators" = "up", 
              "Both" = "both"
            )
          ),
          radioButtons(
            "TF_Filt", "Show only TFs: ", 
            choices = c(
              "Yes" = "TRUE", 
              "No" = "FALSE"
            )
          )
        ),
        
        ##### Main panel ---------------------------------------------------------------------------------------------
        mainPanel(
          plotOutput("sub_GRN", width = "90%", height=400)
        )
      )
    ),
    
    
    ##### TAB -  circuits --------------------------------------------------------------------------------------------
    tabPanel(
      "GRN circuits",
      
      titlePanel("Identifying and extracting 3-node circuits on the basis of an established TF-TF interaction"),
      
      ##### Side bar -------------------------------------------------------------------------------------------
      sidebarLayout(
        sidebarPanel(
          position = "left",width=3,
          selectizeInput('A', 
                         label = "Select upstream TF:",
                         choices = NULL),
          selectizeInput('B', 
                         label = "Select intermediate TF:",
                         choices = NULL)
          
        ),
        
        ##### Main panel ---------------------------------------------------------------------------------------------
        mainPanel(

          fluidRow(
            column(12,h3("FFL:"),DT::dataTableOutput('table1')),
            column(12,h3("Cascade:"),DT::dataTableOutput('table2')),
            column(12,h3("Upstream-TF-specific:"),DT::dataTableOutput('table3'))
          )

        )
      )
    ),
    
    
    
    ##### TAB -  UMAP --------------------------------------------------------------------------------------------
    tabPanel(
      "Patient sub-network clustering",

      titlePanel("Patient sub-network clustering"),
      
      sidebarLayout(
        
        ##### Side bar -------------------------------------------------------------------------------------------
        
        sidebarPanel(
          position = "left", width=3,
          
          radioButtons(
            "BG", "Patient sub-net clustering (UMAP): ", 
            choices = c(
              "Samples" = 1, 
              "Nodes" = 2
            )
          ),
          
          radioButtons(
            "dataset", "Patient sub-net datasets (UMAP):", 
            choices = c(
              "All data" = 1, 
              "ALL" = 2,
              "AML" = 3,
              "FBM" = 4
            )
          ),
          
          downloadButton("downloadData", "Download selection")
          
        ),
        
        ##### Main panel ---------------------------------------------------------------------------------------------
        
        mainPanel(
          fluidRow(column(12), "Use the plotly lasso tool to highlight nodes to download or plot as boxplots (Nodes only)"),
          fluidRow(
            column(6,h3("Dimensionality Reduction:"),plotlyOutput("plot_graph", height=500)),
            column(6,h3("Bin. Expression BoxPlots:"),plotlyOutput("plot_graph2", height=500))
          ),
          fluidRow(
            column(6, 
               radioButtons(
                 "Ann_1", "Sample annotation:", 
                 choices = c(
                   "Gender" = "gender", 
                   "Age" = "age",
                   "Fusion protein" = "FP",
                   "Cytogenic class" = "class",
                   "Dataset" = "group"
                 )
               )
            ),    
            column(6,
               radioButtons(
                 "Ann_2", "Node annotation:", 
                 choices = c(
                   "None" = "none", 
                   "Degree" = "degree",
                   "Stress" = "stress"
                 )
               )
            )
          )
        )
      )
    )
  )
)



