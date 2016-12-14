#' Title
#'
#' @param ExpressionMatrix 
#' @param ModuleList 
#' @param SampleInfo 
#' @param fillMissingValues 
#' @param centerData 
#' @param doubleCenterData 
#' @param minimalNumberOfGenesInModule 
#' @param maximalNumberOfGenesInModule 
#' @param numberOfPermutations 
#' @param typeOfModuleFile 
#'
#' @return
#' @export
#'
#' @examples
rRoma <- function(ExpressionMatrix, ModuleList, SampleInfo=NULL, fillMissingValues = 0, centerData = TRUE, doubleCenterData = FALSE,
                  outlierThreshold = 2, typeOfPCAUsage = 0, robustPCAcalculation = TRUE, robustPCAcalculationForSampling = TRUE, numberOfGeneSetSizesToSample = 5,
                  mostContributingGenesZthreshold = 1, diffSpotGenesZthreshold = 1, correlationThreshold = 0.6, graphicalOutputThreshold = 0.05,
                  minimalNumberOfGenesInModule = 10, maximalNumberOfGenesInModule = 1000, minimalNumberOfGenesInModuleFound = 8,
                  numberOfPermutations = 0, typeOfModuleFile = 0, saveDecomposedFiles = TRUE) {
  
  STANDARD_GMT = 0
  GMT_WITH_WEIGHTS = 1
  WEIGHT_MATRIX = 2
  
  # Getting basic information on the expression matrix
  
  nRow <- nrow(ExpressionMatrix)
  nCol <- ncol(ExpressionMatrix)
  ColNames <- colnames(ExpressionMatrix)
  ColClasses <- rep("", nCol)
  ColDescription <- rep("", nCol)
  RowNames <- unlist(ExpressionMatrix[,1])
  
  ExpressionMatrix <- apply(ExpressionMatrix[,-1], 2, as.numeric)
  nRow <- nrow(ExpressionMatrix)
  nCol <- ncol(ExpressionMatrix)
  
  
  
  print("Constructing basic classes")
  
  TopModuleAnalysis <- new(J("fr.curie.ROMA.ModuleActivityAnalysis"))
  TopModuleAnalysis$table <- new(J("vdaoengine.data.VDataTable"))
  ReadWrite <- .jnew(class = "vdaoengine.data.io.VDatReadWrite")
  TabUtils <- .jnew(class = "fr.curie.ROMA.TableUtils")
  SimpleProcedures <- .jnew(class = "vdaoengine.utils.VSimpleProcedures")
  
  PCAMet <- new(J(class = "vdaoengine.analysis.PCAMethod"))
  PCAMetFC <- new(J(class = "vdaoengine.analysis.PCAMethodFixedCenter"))

  print("Setting up working parameters") 
  
  
  PCAMet$verboseMode <- TRUE
  PCAMetFC$verboseMode <- TRUE
  
  TopModuleAnalysis$fillMissingValues = as.integer(fillMissingValues)
  TopModuleAnalysis$centerData = centerData
  TopModuleAnalysis$doubleCenterData = doubleCenterData
  TopModuleAnalysis$outlierThreshold = .jfloat(outlierThreshold)
  TopModuleAnalysis$typeOfPCAUsage = as.integer(typeOfPCAUsage)
  TopModuleAnalysis$robustPCAcalculation = robustPCAcalculation
  TopModuleAnalysis$robustPCAcalculationForSampling = robustPCAcalculationForSampling
  TopModuleAnalysis$numberOfGeneSetSizesToSample = as.integer(numberOfGeneSetSizesToSample)
  TopModuleAnalysis$mostContributingGenesZthreshold = .jfloat(mostContributingGenesZthreshold)
  TopModuleAnalysis$diffSpotGenesZthreshold = .jfloat(diffSpotGenesZthreshold)
  TopModuleAnalysis$correlationThreshold = .jfloat(correlationThreshold)
  TopModuleAnalysis$graphicalOutputThreshold = .jfloat(graphicalOutputThreshold)
  TopModuleAnalysis$minimalNumberOfGenesInModule = as.integer(minimalNumberOfGenesInModule)
  TopModuleAnalysis$maximalNumberOfGenesInModule = as.integer(maximalNumberOfGenesInModule)
  TopModuleAnalysis$minimalNumberOfGenesInModuleFound = as.integer(minimalNumberOfGenesInModuleFound)
  TopModuleAnalysis$numberOfPermutations = as.integer(numberOfPermutations)
  TopModuleAnalysis$typeOfModuleFile = as.integer(typeOfModuleFile)
  TopModuleAnalysis$saveDecomposedFiles = saveDecomposedFiles
  
  print("Loading gene expression information")
  
  TopModuleAnalysis$table$InjectDoubleTable(as.integer(nRow), as.integer(nCol),
                                            .jarray(ColNames, dispatch = TRUE),
                                            .jarray(ColClasses, dispatch = TRUE),
                                            .jarray(ColDescription, dispatch = TRUE),
                                            .jarray(RowNames, dispatch = TRUE),
                                            .jarray(ExpressionMatrix, dispatch = TRUE))
  
  print("Cleaning up memory")
  rm(ExpressionMatrix)
  gc()
  
  TopModuleAnalysis$table$tableID = as.integer(TopModuleAnalysis$table$lastTableID+1)
  TabUtils$findAllNumericalColumns(TopModuleAnalysis$table)
  
  print("Filling up missing values")
  if(TopModuleAnalysis$fillMissingValues>0){
    print(paste("Filling missing values using first", TopModuleAnalysis$fillMissingValues, "principal components..."))
    TopModuleAnalysis$table = TabUtils$fillMissingValues(TopModuleAnalysis$table, as.integer(fillMissingValues))
  }
  
  
  if(TopModuleAnalysis$centerData){
    print("Centering data (zero mean for any gene)...")
    print("The operation is performed by R")
    
    # Will do in R, the folllwing command will be ignored
    # TopModuleAnalysis$table = TabUtils$normalizeVDat(TopModuleAnalysis$table, TRUE, FALSE)
    
    ExpressionMatrix <- TopModuleAnalysis$table$stringTable
    ColClasses <- TopModuleAnalysis$table$fieldClasses
    ColDescription <- TopModuleAnalysis$table$fieldDescriptions  
    ColNames <- TopModuleAnalysis$table$fieldNames
    RowNames <- ExpressionMatrix[,1]
    
    ExpressionMatrix <- t(scale(t(apply(ExpressionMatrix[,-1], 2, as.numeric)), center = TRUE, scale = FALSE))
    
    TopModuleAnalysis$table$InjectDoubleTable(as.integer(nRow), as.integer(nCol),
                                              .jarray(ColNames, dispatch = TRUE),
                                              .jarray(ColClasses, dispatch = TRUE),
                                              .jarray(ColDescription, dispatch = TRUE),
                                              .jarray(RowNames, dispatch = TRUE),
                                              .jarray(ExpressionMatrix, dispatch = TRUE))
    rm(ExpressionMatrix)
    gc()
    
    TabUtils$findAllNumericalColumns(TopModuleAnalysis$table)
  }
  
  # Double Center Data
  
  if(TopModuleAnalysis$doubleCenterData){
    print("Double centering data (zero mean for any gene and for any sample)...")
    print("The operation is performed by Java. It is probably going to take some time")
    # vd <- .jnew(class = "vdaoengine.data.VDataSet")
    vd <- SimpleProcedures$SimplyPreparedDatasetWithoutNormalization(TopModuleAnalysis$table, as.integer(-1))
    mas <- TabUtils$doubleCenterMatrix(.jarray(vd$massif, dispatch = TRUE))
    
    TopModuleAnalysis$table$InjectDoubleTable(as.integer(nRow), as.integer(nCol),
                                              .jarray(ColNames, dispatch = TRUE),
                                              .jarray(ColClasses, dispatch = TRUE),
                                              .jarray(ColDescription, dispatch = TRUE),
                                              .jarray(RowNames, dispatch = TRUE),
                                              .jarray(apply(.jevalArray(mas, simplify = TRUE), 2, as.double), dispatch = TRUE))
    
    TabUtils$findAllNumericalColumns(TopModuleAnalysis$table)
  }
  
  
  if(TopModuleAnalysis$table$fieldNumByName(TopModuleAnalysis$geneField)==-1){
    TopModuleAnalysis$geneField <- TopModuleAnalysis$table$fieldNames[min(which(TopModuleAnalysis$table$fieldTypes == TopModuleAnalysis$table$STRING))]
  }
  TopModuleAnalysis$table$makePrimaryHash(TopModuleAnalysis$geneField)
  
  ReadGMTDatabase <- function(ModuleList, minimalNumberOfGenesInModule, maximalNumberOfGenesInModule, ReadWeights = TRUE) {
    
    Vec <- new(J("java.util.Vector"))
    
    for(i in 1:length(ModuleList)){
      
      Mg <- new(J("fr.curie.ROMA.Metagene"))
      Mg$name <- ModuleList[[i]]$Name
      Mg$description = ModuleList[[i]]$Desc
      
      GWeight <- rep(NA, length(ModuleList[[i]]$Genes))
      WeiFound <- rep(NA, length(ModuleList[[i]]$Genes))
      
      for(j in 1:length(ModuleList[[i]]$Genes)){
        
        if(length(grep(pattern = "[", x = ModuleList[[i]]$Genes[j], fixed = TRUE))>0){
          GNameWei <- strsplit(ModuleList[[i]]$Genes[j], "[", fixed = TRUE)
          GName = GNameWei[[1]][1]
          GWeight[j] = as.numeric(strsplit(GNameWei[[1]][2], "]", fixed = TRUE)[[1]][1])
          WeiFound[j] <- TRUE
        } else {
          GName <- ModuleList[[i]]$Genes[j]
          GWeight[j] = 1
          WeiFound[j] <- FALSE
        }
        
        Mg$geneNames$add(GName)
        
      }
      
      if(ReadWeights){
        Mg$injectWeights(.jarray(.jfloat(GWeight), dispatch = TRUE), .jarray(WeiFound, dispatch = TRUE))
      }
      
      if((minimalNumberOfGenesInModule<0) | ((Mg$geneNames$size() >= minimalNumberOfGenesInModule) & (Mg$geneNames$size() <= maximalNumberOfGenesInModule))){
        Vec$add(Mg)
      }
    }
    
    return(Vec)
  }
  
  
  print("Processing ModuleList")
  
  sigs = ReadGMTDatabase(ModuleList,minimalNumberOfGenesInModule,maximalNumberOfGenesInModule)
  
  for(i in 1:sigs$size() - 1){
    Mg <- new(J(class = "fr.curie.ROMA.Metagene"), sigs$get(as.integer(i)))
    Mg$probeSets = Mg$geneNames
    
    if(typeOfModuleFile==STANDARD_GMT){
      Mg$initializeWeightsByOnes()
    }
    
    if(typeOfModuleFile==GMT_WITH_WEIGHTS){
      Mg$weights <- sigs$get(as.integer(i))$weights 
      Mg$weightSpecified <- sigs$get(as.integer(i))$weightSpecified 
    }
    
    TopModuleAnalysis$signatures$add(Mg)
  }
  
  
  

  
  ReturnData <- list()
  
  print("Decomposing datasets into modules")
  DataVector <- TopModuleAnalysis$decomposeDataSetDataReturn()
  
  for(i in 1:DataVector$size()){
    ReturnData[[length(ReturnData)+1]] <- .jevalArray(DataVector$get(as.integer(i-1)), simplify = TRUE)
  }
  
  TopModuleAnalysis$numberOfPermutations <- as.integer(numberOfPermutations)
  
  if(TopModuleAnalysis$numberOfPermutations>0){
    print("Calculating sampling")
    explainedVariances_randomDistributions <- TopModuleAnalysis$calcPermutedModuleActivities(TopModuleAnalysis$table)
    DataVector <- TopModuleAnalysis$savePermutedModuleActivitiesDataReturn(explainedVariances_randomDistributions)
    
    for(i in 1:DataVector$size()){
      ReturnData[[length(ReturnData)+1]] <- .jevalArray(DataVector$get(as.integer(i-1)), simplify = TRUE)
    }
  }
  
  print("Compute module activity")
  TopModuleAnalysis$ComputeModuleActivities()
  
  # In the following function I removed a couple of lines!!!
  DataVector <- TopModuleAnalysis$writeGeneProjectionMatrixDataReturn()
  
  for(i in 1:DataVector$size()){
    ReturnData[[length(ReturnData)+1]] <- .jevalArray(DataVector$get(as.integer(i-1)), simplify = TRUE)
  }
  
  print("Produce graphical output")
  DataVector <- TopModuleAnalysis$ProduceGraphicalOutputDataReturn()
  
  for(i in 1:DataVector$size()){
    ReturnData[[length(ReturnData)+1]] <- .jevalArray(DataVector$get(as.integer(i-1)), simplify = TRUE)
  }
  
  print("Top contributing genes determination")
  TopModuleAnalysis$findMostContributingGenes()
  
  print("Creating modile activity table")
  DataVector <- TopModuleAnalysis$moduleActivityTableDataReturn()
  
  for(i in 1:DataVector$size()){
    ReturnData[[length(ReturnData)+1]] <- .jevalArray(DataVector$get(as.integer(i-1)), simplify = TRUE)
  }
  
  print("Cleaning up")
  TopModuleAnalysis$DestroyAll()
  
  for(i in 1:length(ReturnData)){
    names(ReturnData)[i] <- ReturnData[[i]][1,1]
    ReturnData[[i]] <- ReturnData[[i]][-1,]
  }
  
  return(ReturnData)
  
  # sampleFile <- "GTG"
  # 
  # if(!is.null(sampleFile)){
  #   
  #   TopModuleAnalysis$sampleTable  <- new(J("vdaoengine.data.VDataTable"))
  #   
  #   library(readr)
  #   SampleInfo <- read_delim("~/Google Drive/Datasets/Patel et al - Human primary glioblastoma/glioblastoma_samples.txt", 
  #                                      "\t", escape_double = FALSE, trim_ws = TRUE)
  #   
  #   nRow <- nrow(SampleInfo)
  #   nCol <- ncol(SampleInfo)
  #   ColNames <- colnames(SampleInfo)
  #   ColClasses <- rep("", nCol)
  #   ColDescription <- rep("", nCol)
  # 
  #   TopModuleAnalysis$sampleTable$InjectStringVectorAsTable(as.integer(nRow), as.integer(nCol),
  #                                                           .jarray(ColNames, dispatch = TRUE),
  #                                                           .jarray(ColClasses, dispatch = TRUE),
  #                                                           .jarray(ColDescription, dispatch = TRUE),
  #                                                           .jarray(as.character(t(SampleInfo)), dispatch = TRUE))
  #     
  #   TopModuleAnalysis$sampleTable$makePrimaryHash(names(SampleInfo)[1])
  # }
  # 
  # # TopModuleAnalysis$sampleTable$fieldNames
  # 
  # # TopModuleAnalysis$sampleTable$fieldNumByName("ORIGIN")
  # 
  # fieldValueForAveraging <- names(SampleInfo)[3]
  # fieldForAveraging <- names(SampleInfo)[2]
  # 
  # if(is.null(fieldForAveraging)){
  #   System.out.println();
  #   System.out.println("============================================");
  #   System.out.println("     Calculating average module activities");
  #   System.out.println("============================================");
  #   
  #   TopModuleAnalysis$calcAverageModuleActivitiesDataReturn(fieldForAveraging, NULL);
  # }
  # 
  # 
  # TopModuleAnalysis$moduleTable$fieldNames
  # 
  # TopModuleAnalysis$sampleTable$stringTable
  # 
  # 
  # 
  # if(fieldForDiffAnalysis!="NULL"){
  #   if(fieldValuesForDiffAnalysis!=null){
  #     
  #   }
  # }
  #   
  #   
  #   {
  #   System.out.println();
  #   System.out.println("============================================");
  #   System.out.println("     Differential genes determination");
  #   System.out.println("============================================");
  #   
  #   StringTokenizer st = new StringTokenizer(fieldValuesForDiffAnalysis,"#");
  #   String lab1s = st.nextToken();
  #   String lab2s = st.nextToken();
  #   Vector lab1 = new Vector();
  #   Vector lab2 = new Vector();
  #   st = new StringTokenizer(lab1s,",");
  #   while(st.hasMoreTokens())
  #     lab1.add(st.nextToken());
  #   st = new StringTokenizer(lab2s,",");
  #   while(st.hasMoreTokens())
  #     lab2.add(st.nextToken());
  #   findDiffGenes(fieldForDiffAnalysis,lab1,lab2);			
  # }
  
}

