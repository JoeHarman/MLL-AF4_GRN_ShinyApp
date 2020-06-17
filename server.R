rm(list=ls())

###### Load packages ---------------------------------------------------------------------------------------------------------

library(shiny)
library(plotly)
library(tidyverse)
library(rfigshare)
library(igraph)
library(ggnetwork)
library(DT)
library(viridis)
library(intergraph)

shinyServer(function(input, output, session) {

  
  ##### load necessary data ------------------------------------------------------------------------------------------------------
  
  # RNA <- list(
  #   MA4 = read_tsv("RNA-seq_data/contrast_SIMM_vs_SIMA6_cpm.tsv") %>%
  #     mutate(Experiment = "MLL-AF4 siRNA") %>%
  #     select(Experiment, Geneid, logFC, FDR, 
  #            Control_1 = SIMM_1, Control_2 = SIMM_2, Control_3 = SIMM_3, 
  #            Treatment_1 = SIMA6_1, Treatment_2 = SIMA6_2, Treatment_3 = SIMA6_3),
  #   RUNX1 = read_tsv("RNA-seq_data/contrast_NT_NT_vs_RUNX1KD_NT_cpm.tsv") %>%
  #     mutate(Experiment = "RUNX1 siRNA") %>%
  #     select(Experiment, Geneid, logFC, FDR, 
  #            Control_1 = NT_NT_1, Control_2 = NT_NT_2, Control_3 = NT_NT_3, 
  #            Treatment_1 = RUNX1KD_NT_1, Treatment_2 = RUNX1KD_NT_2, Treatment_3 = RUNX1KD_NT_3),
  #   iBET = read_tsv("RNA-seq_data/contrast_DMSO_vs_IBET-1HR_cpm.tsv") %>%
  #     mutate(Experiment = "iBET treatment") %>%
  #     select(Experiment, Geneid, logFC, FDR, 
  #            Control_1 = DMSO.1HR_1, Control_2 = DMSO.1HR_2, Control_3 = DMSO.1HR_3, 
  #            Treatment_1 = IBET.1HR_1, Treatment_2 = IBET.1HR_2, Treatment_3 = IBET.1HR_3),
  #   EPZ = read_tsv("RNA-seq_data/contrast_0umEPZ_vs_2umEPZ_cpm.tsv") %>%
  #     mutate(Experiment = "EPZ treatment") %>%
  #     select(Experiment, Geneid, logFC, FDR, 
  #            Control_1 = X0um.EPZ_1, Control_2 = X0um.EPZ_2, Control_3 = X0um.EPZ_3, 
  #            Treatment_1 = X2um.EPZ_1, Treatment_2 = X2um.EPZ_2, Treatment_3 = X2um.EPZ_3)
  # )
   
  # RNA <- do.call(rbind, RNA) %>%
  #   reshape2::melt(1:4, variable.name = "Samples", value.name = "CPM") %>% 
  #   mutate(Condition = substr(Samples, 1, nchar(as.character(Samples)) - 2)) 
  # 
  # GRN <- list(
  #   edges = read.table("GRN_data/AggregatedGraph_MLL-AF4_edges.tsv", header = T, sep="\t") %>% select(-from, -to),
  #   nodes = read.table("GRN_data/AggregatedGraph_MLL-AF4_nodes.tsv", header = T, sep="\t") %>% select(-name)
  # )
  # DR <- list(
  #   samples = list(
  #     all = read.table("UMAP_data/UMAP_1_binary_all.tsv", header = T, sep="\t"),
  #     ALL = read.table("UMAP_data/UMAP_1_binary_men.tsv", header = T, sep="\t"),
  #     AML = read.table("UMAP_data/UMAP_1_binary_laml.tsv", header = T, sep="\t"),
  #     FBM = read.table("UMAP_data/UMAP_1_binary_fbm.tsv", header = T, sep="\t")
  #   ),
  #   nodes = list(
  #     
  #     all = read.table("UMAP_data/UMAP_2_binary_all.tsv", header = T, sep="\t"),
  #     ALL = read.table("UMAP_data/UMAP_2_binary_men.tsv", header = T, sep="\t"),
  #     AML = read.table("UMAP_data/UMAP_2_binary_laml.tsv", header = T, sep="\t"),
  #     FBM = read.table("UMAP_data/UMAP_2_binary_fbm.tsv", header = T, sep="\t")
  #   )
  # )
   
  # binary.matrix <- data.frame(read.table("UMAP_data//Patient_Binary_Matrix.tsv", header = T, sep="\t"))
  
 # load("bin_mat.RData")
  #load("UMAP.RData")
 # load("RNA.RData")
 # load("GRN.RData")
  
   urltemp <- fs_download(12497765)
   
   
    withProgress(
      message = 'Please wait',
      detail = 'Loading the files...', 
      value = 0, {
       load(url(urltemp[1]))
        incProgress(1/4)
        load(url(urltemp[2]))
        incProgress(2/4)
        load(url(urltemp[3]))
        incProgress(3/4)
        load(url(urltemp[4]))
        incProgress(4/4)
       }
    )

  
  ###### RNA-seq ---------------------------------------------------------------------------------------------------------
  
  updateSelectizeInput(session, 'Gene_Name_RNA', choices = unique(RNA$Geneid), selected="BCL2",server = TRUE)
  
  output$RNA_bp <- renderPlot({
    
    if (input$Gene_Name_RNA=="")
      return("") 
    
    filt_RNA <- RNA %>%
        filter(Geneid == input$Gene_Name_RNA)
    
    ggplot(filt_RNA, aes(x = Condition, y = CPM, col = Condition)) +
      facet_wrap(. ~ Experiment, scales="free") +
      #geom_boxplot() +
      geom_point() +
      expand_limits(y = 0) +
      theme_classic()
    
  })
  

  
  output$table_RNA <- DT::renderDataTable({
    
    if(input$Gene_Name_RNA=="")
      return(NULL)
    
    filt_RNA2 <- RNA %>%
      filter(Geneid == input$Gene_Name_RNA) %>%
      select(Experiment, Geneid, logFC, FDR) %>%
      distinct()

    DT::datatable(filt_RNA2)
  })
  
  
  
  ###### GRN exploratory plots ---------------------------------------------------------------------------------------------------------
  
  updateSelectizeInput(session, "Gene_Name_GRN", choices = unique(GRN$nodes$symbol), selected="KMT2A",server = TRUE)
  
  Filt_GRN <- reactive({
    if (input$Gene_Name_GRN==""){
      return("") 
    }
    
    
    if(input$up_down == "down"){
      e_filt <- GRN$edges %>% filter(from.symbol == input$Gene_Name_GRN)
    }else if(input$up_down == "up"){
      e_filt <- GRN$edges %>% filter(to.symbol == input$Gene_Name_GRN)
    }else{
      e_filt <- GRN$edges %>% filter(to.symbol == input$Gene_Name_GRN | from.symbol == input$Gene_Name_GRN)
    }
    
    g_filt <- unique(c(as.character(e_filt$to.symbol), as.character(e_filt$from.symbol)))
    
    n_filt <- GRN$nodes %>% 
      filter(symbol %in% g_filt) %>%
      distinct(symbol, .keep_all = T)
    
    if(input$TF_Filt == "TRUE"){
      n_filt <- n_filt %>% filter(stress > 1)
      e_filt <- e_filt %>% filter(to.symbol %in% n_filt$symbol & from.symbol %in% n_filt$symbol)
    }
    
    return(list(e_filt, n_filt))
    
  })
  
  output$sub_GRN <- renderPlot({
    
    
    if(input$Gene_Name_GRN=="")
      return("") 
    
    
    Filt_GRN2 <- Filt_GRN()
    
    #if(nrow(Filt_GRN2[[1]]) > 500)
    #  return("Too many interactions")
    
    g <- igraph::graph_from_data_frame(Filt_GRN2[[1]], T, Filt_GRN2[[2]])
    fortified.graph <- g %>% 
      ggnetwork::ggnetwork(layout = "kamadakawai")
    
    fortified.graph %>% 
      ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(
        color = "grey60", 
        size = .15, curvature = 0.1, 
        arrow = arrow(length = unit(6, "pt"), type = "open")
      ) +
      geom_nodes(
        aes(fill = degree),
        pch=21,
        size=5
      ) +
      geom_nodetext(aes(label = vertex.names), col="black") +
      scale_fill_viridis() +
      theme_blank() +
      theme(legend.position = "bottom")
    
    
    
  })
  
  ###### GRN circuit plots ---------------------------------------------------------------------------------------------------------
  
  updateSelectizeInput(session, "A", choices = unique(GRN$nodes$symbol[GRN$nodes$stress > 1]), selected="KMT2A",server = TRUE)
  updateSelectizeInput(session, "B", choices = unique(GRN$nodes$symbol[GRN$nodes$stress > 1]), selected="RUNX1",server = TRUE)
  
  Filt_circuit <- reactive({
    
    if (input$A==""){
      return("") 
    }
    if (input$B==""){
      return("") 
    }
    
    tmp <- paste(GRN$edges$from.symbol, GRN$edges$to.symbol, sep="_")
    
    if(!(paste(input$A, input$B, sep="_") %in% tmp)){
      return("Not an existing A-B interaction")
    }
    
    e_A_target <- GRN$edges$to.symbol[GRN$edges$from.symbol == input$A & !(GRN$edges$to.symbol %in% c(input$A, input$B)) ]
    e_B_target <- GRN$edges$to.symbol[GRN$edges$from.symbol == input$B & !(GRN$edges$to.symbol %in% c(input$A, input$B)) ]
    
    g.ffl <- intersect(e_A_target, e_B_target)
    g.cascade <- setdiff(e_B_target, e_A_target)
    g.A_specific <- setdiff(e_A_target, e_B_target)
    
    n.ffl <- GRN$nodes %>%
      filter(symbol %in% g.ffl) %>%
      distinct(symbol, .keep_all = T) %>%
      data.frame() %>%
      select(symbol, degree, stress, Dropout_class)
    n.cascade <- GRN$nodes %>%
      filter(symbol %in% g.cascade) %>%
      distinct(symbol, .keep_all = T) %>%
      data.frame() %>%
      select(symbol, degree, stress, Dropout_class)
    n.A_specific <- GRN$nodes %>%
      filter(symbol %in% g.A_specific) %>%
      distinct(symbol, .keep_all = T) %>%
      data.frame() %>%
      select(symbol, degree, stress, Dropout_class)
    
  
    return(list(n.ffl, n.cascade, n.A_specific))
    
  })
  

  
  output$table1 <- DT::renderDataTable(DT::datatable(Filt_circuit()[[1]]))
  output$table2 <- DT::renderDataTable(DT::datatable(Filt_circuit()[[2]]))
  output$table3 <- DT::renderDataTable(DT::datatable(Filt_circuit()[[3]]))
  
  ###### UMAP plots ---------------------------------------------------------------------------------------------------------
  
  ### Select UMAP dataset ###
  x <- reactive({
    DR[[as.numeric(input$BG)]][[as.numeric(input$dataset)]]
  })
  
  
  ### Plot UMAP ###
  output$plot_graph <- renderPlotly({
    
    a <- 0.8 # Alpha for points
    
    # If clustering samples...
    if(input$BG == "1"){
      p <- ggplot(x(), aes(x = PC1, y = PC2, text= paste0("ID: ", sample, "\n Gender: ", gender,"\n Age: ", age,"\n FP: ", FP,"\n class: ", class)))
        if(input$Ann_1 == "gender"){p <- p + geom_point(aes(col = gender), alpha=a)}
        if(input$Ann_1 == "age"){p <- p + geom_point(aes(col = age), alpha=a)}
        if(input$Ann_1 == "FP"){p <- p + geom_point(aes(col = FP), alpha=a)}
        if(input$Ann_1 == "class"){p <- p + geom_point(aes(col = class), alpha=a)}
        if(input$Ann_1 == "group"){p <- p + geom_point(aes(col = group), alpha=a)}
        p <- p + theme_classic()
      
    # Else, if clustering nodes...
    }else{
      p <- ggplot(x()[order(x()$degree),], aes(x = PC1, y = PC2, text=paste0("Node: ", sample, "\n Stress: ", stress)))
        if(input$Ann_2 == "none"){p <- p + geom_point(col = "lightblue", alpha=a)}
        if(input$Ann_2 == "degree"){p <- p + geom_point(aes(col = degree), alpha=a) + scale_color_gradient(low="grey70", high="red")}
        if(input$Ann_2 == "stress"){p <- p + geom_point(aes(col = stress), alpha=a) + scale_color_gradient(low="grey70", high="red")}
        p <-  p + theme_classic()
    }
    
    # Create plotly object
    plotme <- ggplotly(p, originalData=T)
    
    # Register selections in plotly
    event_register(plotme, "plotly_selected")
    
  })
  
  
  ### Plot UMAP selection boxplots ###
  output$plot_graph2 <- renderPlotly({
    
    # Use registered plotly selection
    d <- event_data("plotly_selected")

    df <- x()
    df <- df$sample[signif(df$PC1, 4) %in% signif(d$x, 4) & signif(df$PC2, 4) %in% signif(d$y, 4)]
    
    # If clustering samples...
    if(input$BG == "1"){
      
      # No plot
      p2 <- ggplot()
      
    # If clustering nodes...
    }else{
      
      # Split binary matrix by dataset
      binary.matrix.men <- binary.matrix[binary.matrix[,1] %in% df,c(1, 375:429)]
      binary.matrix.laml <- binary.matrix[binary.matrix[,1] %in% df,c(1, 196:374)]
      binary.matrix.ox_fbm <- binary.matrix[binary.matrix[,1] %in% df,c(1, 430:460)]
      
      # Long format
      binary.matrix.long <- rbind(
        data.frame(Geneid = binary.matrix.men[,1], mean_binary = apply(binary.matrix.men[,-1], 1, mean), group="ALL"),
        data.frame(Geneid = binary.matrix.laml[,1], mean_binary = apply(binary.matrix.laml[,-1], 1, mean), group="AML"),
        data.frame(Geneid = binary.matrix.ox_fbm[,1], mean_binary = apply(binary.matrix.ox_fbm[,-1], 1, mean), group="FBM")
        
      )
      
      # Create boxplots based on data
      p2 <- ggplot(binary.matrix.long, aes(x = group, y = mean_binary, fill = group)) +
        geom_boxplot() +
        expand_limits(x = 0, y = 0) +
        theme_classic()
    
      
    }
    
    if(input$BG == "2"){ggplotly(p2)}else{(p2)}
  })
  
  
  ### Download plotly selection ###
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste("Download", ".csv", sep = "")
    },
    
    content = function(file) {
      d <- event_data("plotly_selected")
      df <- x()
      df <- df[signif(df$PC1, 4) %in% signif(d$x, 4) & signif(df$PC2, 4) %in% signif(d$y, 4),]
      
      write.csv(df, file, row.names = FALSE)
    }
    
  )
    
})
