#' Title
#'
#' @param RomaData 
#' @param Groups 
#' @param ExpMat 
#'
#' @return
#' @export
#'
#' @examples
DashboardA <- function(RomaData, Groups, ExpMat, Interactive = TRUE) {
  
  #
  # This is a Shiny web application. You can run the application by clicking
  # the 'Run App' button above.
  #
  # Find out more about building applications with Shiny here:
  #
  #    http://shiny.rstudio.com/
  #
  
  library(shiny)
  if(Interactive){
    library(plotly)
  } else {
    if(R.utils::isPackageLoaded("plotly")){
      detach("package:plotly", unload=TRUE)
    }
  }
 
  ProcessedSamples <- colnames(RomaData$ProjMatrix)
  
  if(is.null(Groups)){
    Groups <- rep("N/A", length(ProcessedSamples))
    names(Groups) <- ProcessedSamples
  }
  
  GSList <- as.list(1:nrow(RomaData$PVVectMat))
  names(GSList) <- as.list(rownames(RomaData$ProjMatrix))
  GSList <- GSList[order(names(GSList))]
  
  PCAProj <- irlba::prcomp_irlba(t(ExpMat[,ProcessedSamples]), 2, retx = TRUE)$x 
  rownames(PCAProj) <- ProcessedSamples
  colnames(PCAProj) <- c("PC1", "PC2")
  
  # Projs <- PCAProj
  
  # Define UI for application that draws a histogram
  ui <- pageWithSidebar(
    
    # Application title
    headerPanel("rRoma / Sample dashboard"),
    
    # Sidebar with a slider input for number of bins
    
    sidebarPanel(
      selectInput("prjt", "Projectin type:",
                  list("PCA" = "PCA", "tSNE" = "tSNE")),
      sliderInput("perp", "tSNE perplexity:",
                  min = 0,  max = 50,  value = 0, step = .1),
      hr(),
      
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
      if(Interactive){
        plotlyOutput("scatPlot")
      } else {
        plotOutput("scatPlot")
      }
      
    )
    
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
    
    SelectedGS <- reactive({
      as.integer(input$gs)
    })
    
    GetProj <- reactive({
      
      if(input$prjt == "PCA"){
        PCAProj
      }
      
      if(input$prjt == "tSNE"){
        tSNEProj <- Rtsne::Rtsne(X = t(ExpMat[,ProcessedSamples]), perplexity = input$perp)$Y
        rownames(tSNEProj) <- ProcessedSamples
        colnames(tSNEProj) <- c("PC1", "PC2")
        
        tSNEProj
      }
    })
    
    observe({
      
      Idx <- SelectGeneSets(RomaData = RomaData, VarThr = 10^input$pdisp, VarMode = "Wil",
                            VarType = input$disp,
                            MedThr = 10^input$pexp, MedMode = "Wil",
                            MedType = input$exp)
      
      GSList <- Idx
      names(GSList) <- as.list(rownames(RomaData$ProjMatrix)[Idx])
      GSList <- GSList[order(names(GSList))]
      
      updateSelectInput(session, "gs", choices = GSList, selected = GSList[[1]])
      
    })
    
    output$info <- renderUI({
      str1 <- paste("Overdispersed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 1]))
      str2 <- paste("Underdispersed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 2]))
      str3 <- paste("Overexpressed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 5]))
      str4 <- paste("Underexpressed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 6]))
      HTML(paste(str1, str2, str3, str4, sep = '<br/>'))
      
    })
    
    if(Interactive){
      output$scatPlot <- renderPlotly({
        
        Projs <- GetProj()
        
        if(is.null(Projs)){
          Projs <- PCAProj
        }
        
        p <- ggplot2::ggplot(data = data.frame(Comp1 = Projs[,"PC1"],
                                               Comp2 = Projs[,"PC2"],
                                               Score = RomaData$ProjMatrix[SelectedGS(), ProcessedSamples],
                                               Group = Groups[ProcessedSamples]),
                             ggplot2::aes(x = Comp1, y = Comp2, shape = Group, color = Score)) +
          ggplot2::geom_point(ggplot2::aes(text = rownames(Projs))) +
          ggplot2::scale_color_gradient2(low = "blue", high = "red", mid = "white") +
          ggplot2::labs(x = "Component 1", y = "Component 2", shape = "",
                        title = rownames(RomaData$ProjMatrix)[SelectedGS()])
        
        print(ggplotly(p))
      })
    } else {
      output$scatPlot <- renderPlot({
        
        Projs <- GetProj()
        
        if(is.null(Projs)){
          Projs <- PCAProj
        }
        
        p <- ggplot2::ggplot(data = data.frame(Comp1 = Projs[,"PC1"],
                                               Comp2 = Projs[,"PC2"],
                                               Score = RomaData$ProjMatrix[SelectedGS(), ProcessedSamples],
                                               Group = Groups[ProcessedSamples]),
                             ggplot2::aes(x = Comp1, y = Comp2, shape = Group, color = Score)) +
          ggplot2::geom_point(ggplot2::aes(text = rownames(Projs)), size = 3) +
          ggplot2::scale_color_gradient2(low = "blue", high = "red", mid = "white") +
          ggplot2::labs(x = "Component 1", y = "Component 2", shape = "",
                        title = rownames(RomaData$ProjMatrix)[SelectedGS()])
        
        print(p)
      })
    }
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
  
}
