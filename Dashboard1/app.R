DashboardA <- function(RomaData, Groups, ExpMat) {
  
  #
  # This is a Shiny web application. You can run the application by clicking
  # the 'Run App' button above.
  #
  # Find out more about building applications with Shiny here:
  #
  #    http://shiny.rstudio.com/
  #
  
  library(shiny)
  
  ProcessedSamples <- colnames(RomaData$ProjMatrix)
  
  GSList <- as.list(1:nrow(RomaData$PVVectMat))
  names(GSList) <- as.list(rownames(RomaData$ProjMatrix))
  GSList <- GSList[order(names(GSList))]
  
  Proj <- irlba::prcomp_irlba(t(ExpMat[,ProcessedSamples]), 2, retx = TRUE)$x 
  rownames(Proj) <- ProcessedSamples
  colnames(Proj)
  
  # Define UI for application that draws a histogram
  ui <- pageWithSidebar(
    
    # Application title
    headerPanel("rRoma dashboard A"),
    
    # Sidebar with a slider input for number of bins
    
    sidebarPanel(
      selectInput("disp", "Dispersion filter:",
                  list("Overdispersed" = "Over", "Underdispersed" = "Under", "None" = "None")),
      sliderInput("pdisp", "Log 10 p-value threshold:",
                  max = 0,  min = -10,  value = 0, step = .1),
      hr(),
      
      selectInput("exp", "Expression filter:",
                  list("Overexpressed" = "Over", "Underexpressed" = "Under", "None" = "None")),
      sliderInput("pexp", "Log 10 p-value threshold:",
                  max = 0,  min = -10,  value = 0, step = .1),
      hr(),
      
      selectInput("gs", "GeneSet:",
                  GSList),
      htmlOutput("info")
    ),
    
    mainPanel(
      plotly::plotlyOutput("scatPlot")
    )
    
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
    
    SelectedGS <- reactive({
      as.integer(input$gs)
    })
    
    
    observe({
      
      # print(input$disp)
      # print(input$exp)
      
      
      Idx <- SelectGeneSets(RomaData = RomaData, VarThr = 10^input$pdisp, VarMode = "Wil",
                            VarType = input$disp,
                            MedThr = 10^input$pexp, MedMode = "Wil",
                            MedType = input$exp)
      
      # print(Idx)
      
      # SelectGeneSets(RomaData = RomaData, VarThr = 1, VarMode = "Wil",
      #                VarType = "Over",
      #                MedThr = 1, MedMode = "Wil",
      #                MedType = "input$exp")
      
      GSList <- Idx
      names(GSList) <- as.list(rownames(RomaData$ProjMatrix)[Idx])
      GSList <- GSList[order(names(GSList))]
      
      updateSelectInput(session, "gs", choices = GSList, selected = Idx[1])
      
      # choice1<-c("NONE",setdiff(heirarchy,c(hei2,hei3)))
      # choice2<-c("NONE",setdiff(heirarchy,c(hei1,hei3)))
      # choice3<-c("NONE",setdiff(heirarchy,c(hei1,hei2)))
      # 
      # updateSelectInput(session,"heir1",choices=choice1,selected=hei1)
      # updateSelectInput(session,"heir2",choices=choice2,selected=hei2)
      # updateSelectInput(session,"heir3",choices=choice3,selected=hei3)
      
    })
    
    
    output$info <- renderUI({
      str1 <- paste("Overdispersed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 1]))
      str2 <- paste("Underdispersed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 2]))
      str3 <- paste("Overexpressed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 5]))
      str4 <- paste("Underexpressed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 6]))
      HTML(paste(str1, str2, str3, str4, sep = '<br/>'))
      
    })
    
    output$scatPlot <- plotly::renderPlotly({
      
      p <- ggplot2::ggplot(data = data.frame(Comp1 = Proj[,"PC1"],
                                             Comp2 = Proj[,"PC2"],
                                             Score = RomaData$ProjMatrix[SelectedGS(), ProcessedSamples],
                                             Group = Groups[ProcessedSamples]),
                           ggplot2::aes(x = Comp1, y = Comp2, shape = Group, color = Score)) +
        ggplot2::geom_point(ggplot2::aes(text = rownames(Proj))) +
        ggplot2::scale_color_gradient2(low = "blue", high = "red", mid = "white") +
        ggplot2::labs(x = "Component 1", y = "Component 2", shape = "")
      
      print(plotly::ggplotly(p))
      # print(p)
    })
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
  
}
