options("java.home"="/Library/Java/JavaVirtualMachines/jdk1.8.0_102.jdk/Contents/Home/jre")
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_102.jdk/Contents/Home/jre/lib/server/libjvm.dylib')

options(java.parameters = "-Xmx4g")

library(rRoma)

# Load Data

library(readr)

ExpressionMatrix <- read_delim("~/Google Drive/Datasets/Patel et al - Human primary glioblastoma/PATEL_uniquelog.zip", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

ExpressionMatrix <- ExpressionMatrix[, -ncol(ExpressionMatrix)]

GeneNames <- unlist(ExpressionMatrix[,1])
# ExpressionMatrix <- ExpressionMatrix[,-ncol(ExpressionMatrix)]
# colnames(ExpressionMatrix) <- "GeneNames"

# library(GO.db)

# library(rpgraph)

ModuleList <- list(list(Name="G1/S", Desc="G1/S",
                        Genes = intersect(GeneNames, union(rpgraph::StageAssociation_Whit$S1_U, rpgraph::StageAssociation_Whit$S2_U))),
                   list(Name="G2/M", Desc="G2/M",
                        Genes = intersect(GeneNames, union(rpgraph::StageAssociation_Whit$S3_U, rpgraph::StageAssociation_Whit$S4_U)))
)

SampleInfo <- read_delim("~/Google Drive/Datasets/Patel et al - Human primary glioblastoma/glioblastoma_samples.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

ToKeep <- SampleInfo$SAMPLE[SampleInfo$TYPE == "Single cell mRNA-seq"]

ToKeepID = c(1, which(colnames(ExpressionMatrix) %in% ToKeep))
ExpressionMatrix <- ExpressionMatrix[,ToKeepID]

Tables <- rRoma(ExpressionMatrix = ExpressionMatrix, ModuleList = ModuleList, numberOfPermutations = 50)

plot(t(Tables$moduletable_simple.txt[-1, -1]),
     xlab = Tables$moduletable_simple.txt[2,1],
     ylab = Tables$moduletable_simple.txt[3,1])
