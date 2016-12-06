library(rJava)

.jinit()
.jaddClassPath(path = "/bioinfo/guests/lalberga/R/rRoma/Roma-master/jar/ROMA.jar")
# Load Data

ExampleData <- matrix(as.character(abs(rnorm(10*100))), 10)

# Constructing VDataTable

VDataTable <- .jnew(class = "vdaoengine/data/VDataTable")

VDataTable$stringTable <- .jarray(ExampleData, dispatch = TRUE)

, new.class = "[[Ljava/lang/String;")
