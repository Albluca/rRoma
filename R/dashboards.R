#' Sample dashboard
#'
#' @param RomaData 
#' @param Groups 
#' @param ExpMat 
#' @param Interactive 
#'
#' @return
#' @export
#'
#' @examples
#' 
DashboardA <- function(RomaData, ExpMat, Groups = NULL, Interactive = FALSE) {
  
  library(shiny)
  
  if(Interactive){
    print("Using plotly. This can cause problems on some systems. Try setting 'Interactive = FALSE' if errors are encountered")
    library(plotly)
  } else {
    if(R.utils::isPackageLoaded("plotly")){
      print("Detaching plotly.")
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
        tabPanel("Plot",
                 plotlyOutput("scatPlot"),
                 plotOutput("boxPlot")
        )
      } else {
        tabPanel("Plot",
                 plotOutput("scatPlot"),
                 plotOutput("boxPlot")
        )
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
    
    output$boxPlot <- renderPlot({
      
      Projs <- GetProj()
      
      if(is.null(Projs)){
        Projs <- PCAProj
      }
      
      GetComb <- function(GrpLevs) {
        RetList <- list()
        for(i in 1:length(GrpLevs)){
          for(j in 1:length(GrpLevs)){
            if(i<j){
              RetList[[length(RetList)+1]] <- c(i, j)
            }
          }
        }
        return(RetList)
      }
      
      p <- ggplot2::ggplot(data = data.frame(Score = RomaData$ProjMatrix[SelectedGS(), ProcessedSamples],
                                             Group = Groups[ProcessedSamples]),
                           ggplot2::aes(x = Group, y = Score)) +
        ggplot2::geom_boxplot() +
        ggsignif::geom_signif(comparisons = GetComb(unique(Groups[ProcessedSamples])),
                              map_signif_level=TRUE, test = "wilcox.test", step_increase = .1)
      
      print(p)
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
          ggplot2::geom_point(size = 3) +
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
























#' Genesets dashboard
#'
#' @param RomaData 
#' @param Groups 
#'
#' @return
#' @export
#'
#' @examples
#' 
DashboardB <- function(RomaData, Groups = NULL) {
  
  library(shiny)
  
  ProcessedSamples <- colnames(RomaData$ProjMatrix)
  
  if(!is.null(Groups)){
    
    FoundSamp <- intersect(colnames(RomaData$ProjMatrix), names(Groups))
    Groups <- Groups[FoundSamp]
    
    if(length(FoundSamp) > 2){
      if(!is.factor(Groups)){
        Groups <- as.factor(Groups)
      }
      AddInfo <- data.frame(Groups = Groups)
    } else {
      AddInfo = NULL
    }
    SelList <- list("By group" = "group", "By sample" = "sample")
    SelListAF <- list("Mean" = "mean", "Median" = "median", "Std. dev." = "sd", "IQR" = "IQR")
    
  } else {
    AddInfo = NULL
    SelList <- list("By sample" = "sample")
    SelListAF <- list()
  }
  
  
  # Projs <- PCAProj
  
  # Define UI for application that draws a histogram
  ui <- pageWithSidebar(
    
    # Application title
    headerPanel("rRoma / Genesets dashboard"),
    
    # Sidebar with a slider input for number of bins
    
    sidebarPanel(
      
      selectInput("htype", "Heatmap type:", SelList),
      selectInput("aggfun", "Aggregating function:", SelListAF),
      hr(),
      
      checkboxInput("gsclus", "Cluster genesets", FALSE),
      checkboxInput("saclus", "Cluster samples / groups", FALSE),
      checkboxInput("gscol", "Samples on columns", FALSE),
      hr(),
      
      selectInput("disp", "Dispersion filter:",
                  list("Overdispersed" = "Over", "Underdispersed" = "Under", "None" = "None")),
      sliderInput("pdisp", "Log 10 p-value threshold:",
                  max = 0,  min = -10,  value = 0, step = .1),
      hr(),
      
      selectInput("exp", "Expression filter:",
                  list("Overexpressed" = "Over", "Underexpressed" = "Under", "None" = "None")),
      sliderInput("pexp", "Log 10 p-value threshold:",
                  max = 0,  min = -10,  value = 0, step = .1)
      
    ),
    
    mainPanel(
      plotOutput("hmPlot", height = "800px")
    )
    
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {

    SelGS <- reactive({
      SelectGeneSets(RomaData = RomaData, VarThr = 10^input$pdisp, VarMode = "Wil",
                     VarType = input$disp,
                     MedThr = 10^input$pexp, MedMode = "Wil",
                     MedType = input$exp)
    })
    
    ClusterGS <- reactive({
      input$gsclus
    })
    
    ClusterSA <- reactive({
      input$saclus
    })
    
    Transpose <- reactive({
      input$gscol
    })
    
    
    output$hmPlot <- renderPlot({
      
      SelectedGS <- SelGS()
      
      PlotMat <- RomaData$ProjMatrix[SelectedGS,]
      
      if(input$htype == "sample"){
        
        if(length(SelectedGS)>1){
          
          if(Transpose()){
            pheatmap::pheatmap(PlotMat, color = colorRamps::blue2red(50),
                               cluster_rows = ClusterGS(), cluster_cols = ClusterSA(),
                               annotation_col = AddInfo)
          } else {
            pheatmap::pheatmap(t(PlotMat), color = colorRamps::blue2red(50),
                               cluster_rows = ClusterSA(), cluster_cols = ClusterGS(),
                               annotation_row = AddInfo)
          }
          
          
        }
        
        if(length(SelectedGS) == 1){
          
          names(PlotMat) <- colnames(RomaData$ProjMatrix)
          
          pheatmap::pheatmap(t(PlotMat), color = colorRamps::blue2red(50),
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             main = rownames(RomaData$ProjMatrix)[SelectedGS])
          
        }
      }
      
      if(input$htype == "group"){
        
        if(length(SelectedGS) > 1){
          
          SplitData <- split(data.frame(t(PlotMat[,FoundSamp])), f=AddInfo$Groups)
          
          Aggmat <- sapply(SplitData, function(x) {
            apply(x, 2, get(input$aggfun))
          })
          
          if(Transpose()){
            pheatmap::pheatmap(Aggmat, color = colorRamps::blue2red(50),
                               cluster_rows = ClusterGS(), cluster_cols = ClusterSA())
          } else {
            pheatmap::pheatmap(t(Aggmat), color = colorRamps::blue2red(50),
                               cluster_rows = ClusterSA(), cluster_cols = ClusterGS())
          }
          
        }
        
        if(length(SelectedGS) == 1){
          
          # names(PlotMat) <- colnames(RomaData$ProjMatrix)
          SplitData <- split(data.frame(PlotMat[FoundSamp]), f=AddInfo$Groups)
          
          Aggmat <- sapply(SplitData, function(x) {
            do.call(input$aggfun, list(unlist(x)))
          })
          
          pheatmap::pheatmap(t(Aggmat), color = colorRamps::blue2red(50),
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             main = rownames(RomaData$ProjMatrix)[SelectedGS])
        }
        
      }
      
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
  
}
