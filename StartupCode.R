library(rJava)

.jinit()
.jaddClassPath(path = "/bioinfo/guests/lalberga/R/rRoma/Roma-master/jar/ROMA.jar")
# Load Data

ExampleData <- matrix(as.character(abs(rnorm(10*100))), 10)

# Constructing VDataTable

Tabl <- .jnew(class = "vdaoengine/data/VDataTable")
Tabl$colCount <- ncol(ExampleData)
Tabl$fieldNames <- .jarray(paste("N", 1:ncol(ExampleData), sep = ''), "java/lang/String", dispatch = TRUE)


VDataTable$stringTable <- .jarray(ExampleData, dispatch = TRUE)

, new.class = "[[Ljava/lang/String;")
