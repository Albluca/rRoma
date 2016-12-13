# install.packages(c("roxygen2", "testthat", "knitr", "rstudioapi"))
# devtools::install_github("hadley/devtools")

# install.packages(c("formatR", "lintr"))

#
# devtools::use_package("rJava", "Depends")
#
# devtools::use_package("igraph", "Imports")
# devtools::use_package("tictoc", "Imports")
# devtools::use_package("rgl", "Imports")
#
# devtools::use_package("bigpca", "Suggests")
# devtools::use_package("flashpcaR", "Suggests")
# devtools::use_package("irlba", "Suggests")
# devtools::use_package("nsprcomp", "Suggests")
# devtools::use_package("repmis", "Suggests")
# devtools::use_package("plotly", "Suggests")
# devtools::use_package("devtools", "Suggests")
# devtools::use_package("repmis", "Suggests")
# devtools::use_package("downloader", "Suggests")
#
# devtools::use_package("pcaMethods", "Suggests")
#
#
# 
# 
# 
# 
# load("~/Datasets/Human_Gene_Info.RData")
# 
# library(GO.db)
# AllTerms <- GOBPOFFSPRING[["GO:0007049"]]
# HumanGenes_GOCellCycle <- unique(Human_Gene_Info$hgnc_symbol[which(Human_Gene_Info$go_id %in% AllTerms)])
# 
# 
# AllTerms <- GOBPOFFSPRING[["GO:0006260"]]
# HumanGenes_GODNAReplication <- unique(Human_Gene_Info$hgnc_symbol[which(Human_Gene_Info$go_id %in% AllTerms)])
# 
# 
# load("~/Datasets/Mouse_Gene_Info.RData")
# 
# library(GO.db)
# AllTerms <- GOBPOFFSPRING[["GO:0007049"]]
# MouseGenes_GOCellCycle <- unique(Mouse_Gene_Info$external_gene_name[which(Mouse_Gene_Info$go_id %in% AllTerms)])
# 
# 
# AllTerms <- GOBPOFFSPRING[["GO:0006260"]]
# MouseGenes_GODNAReplication <- unique(Mouse_Gene_Info$external_gene_name[which(Mouse_Gene_Info$go_id %in% AllTerms)])
# 
# 
# 
# 
# library(devtools)
# use_data(HumanGenes_GOCellCycle)
# use_data(MouseGenes_GOCellCycle)
# use_data(HumanGenes_GODNAReplication)
# use_data(MouseGenes_GODNAReplication)
# 
# 
# 
# 
# 
# library(downloader)
#
# UrlLoc <- "https://raw.githubusercontent.com/auranic/Elastic-principal-graphs/master/ElasticPrincipalGraphs_matlab/test_data/circle/simple_circle.data"
# TempLoc <- tempfile(pattern = "simple_circle", fileext = ".data")
# download.file(url = UrlLoc, destfile = TempLoc, mode = "w")
# simple_circle <- read.delim(TempLoc, header=FALSE, stringsAsFactors=FALSE)
# file.remove(TempLoc)
#
# devtools::use_data(simple_circle)
#
#
#
# UrlLoc <- "https://raw.githubusercontent.com/auranic/Elastic-principal-graphs/master/ElasticPrincipalGraphs_matlab/test_data/tree23/tree23.data"
# TempLoc <- tempfile(pattern = "tree23", fileext = ".data")
# download.file(url = UrlLoc, destfile = TempLoc, mode = "w")
# simple_tree <- read.delim(TempLoc, header=FALSE, stringsAsFactors=FALSE)
# file.remove(TempLoc)
#
# devtools::use_data(simple_tree)
#
#
#
# UrlLoc <- "https://raw.githubusercontent.com/auranic/Elastic-principal-graphs/master/ElasticPrincipalGraphs_matlab/test_data/iris/iris.txt"
# TempLoc <- tempfile(pattern = "iris", fileext = ".txt")
# download.file(url = UrlLoc, destfile = TempLoc, mode = "w")
# simple_iris <- read.delim(TempLoc, header=TRUE, stringsAsFactors=FALSE)
# file.remove(TempLoc)
#
# devtools::use_data(simple_iris)

# devtools::use_vignette("Simple_Examples")









# 
# 
# 
# WHITFIELD_CELL_CYCLE_S <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_S.txt")
# WHITFIELD_CELL_CYCLE_G1_S <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_G1_S.txt")
# WHITFIELD_CELL_CYCLE_G2_M <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_G2_M.txt")
# WHITFIELD_CELL_CYCLE_M_G1 <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_M_G1.txt")
# WHITFIELD_CELL_CYCLE_G2 <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_G2.txt")
# 
# 
# 
# StageAssociation_Whit <- list(QVarCutOff = 0.7,
#                               Stages = c("G1", "S", "G2", "M"),
#                               S1_U = unique(c(WHITFIELD_CELL_CYCLE_G1_S[[1]],
#                                               WHITFIELD_CELL_CYCLE_M_G1[[1]])),
#                               S2_U = unique(c(WHITFIELD_CELL_CYCLE_G1_S[[1]],
#                                               WHITFIELD_CELL_CYCLE_S[[1]])),
#                               S3_U = unique(c(WHITFIELD_CELL_CYCLE_G2[[1]],
#                                               WHITFIELD_CELL_CYCLE_G2_M[[1]])),
#                               S4_U = unique(c(WHITFIELD_CELL_CYCLE_G2_M[[1]],
#                                               WHITFIELD_CELL_CYCLE_M_G1[[1]])))
# 
# 
# devtools::use_data(StageAssociation_Whit)
# 
# StageAssociation_Whit_G0 <- StageAssociation_Whit
# StageAssociation_Whit_G0$Stages <- c("G0", StageAssociation_Whit$Stages)
# StageAssociation_Whit_G0$S1_D <- unique(c(StageAssociation_Whit$S1_U, StageAssociation$S2_U, StageAssociation$S3_U, StageAssociation$S4_U))
# StageAssociation_Whit_G0$S1_U <- NULL
# StageAssociation_Whit_G0$S2_U <- StageAssociation_Whit$S1_U
# StageAssociation_Whit_G0$S3_U <- StageAssociation_Whit$S2_U
# StageAssociation_Whit_G0$S4_U <- StageAssociation_Whit$S3_U
# StageAssociation_Whit_G0$S5_U <- StageAssociation_Whit$S4_U
# 
# 
# devtools::use_data(StageAssociation_Whit_G0)
# 
# 
# library(readr)
# 
# WHITFIELD_CELL_CYCLE_S <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_S.txt")
# WHITFIELD_CELL_CYCLE_G1_S <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_G1_S.txt")
# WHITFIELD_CELL_CYCLE_G2_M <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_G2_M.txt")
# WHITFIELD_CELL_CYCLE_M_G1 <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_M_G1.txt")
# WHITFIELD_CELL_CYCLE_G2 <- read_csv("/bioinfo/users/lalberga/Datasets/FromAndrei/WHITFIELD_CELL_CYCLE_G2.txt")
# 
# 
# StageAssociation_Whit_Ext <- list(QVarCutOff = 0.7,
#                               Stages = c("G1/S", "S", "G2", "G2/M", "M/G1"),
#                               S1_U = WHITFIELD_CELL_CYCLE_G1_S[[1]],
#                               S2_U = WHITFIELD_CELL_CYCLE_S[[1]],
#                               S3_U = WHITFIELD_CELL_CYCLE_G2[[1]],
#                               S4_U = WHITFIELD_CELL_CYCLE_G2_M[[1]],
#                               S5_U = WHITFIELD_CELL_CYCLE_M_G1[[1]])
# 
# 
# devtools::use_data(StageAssociation_Whit_Ext)
# 
