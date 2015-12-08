

list.of.packages <- c(

"biomaRt",
"reshape2",
"minfi",
"rtracklayer",
"rmarkdown",
"markdown",
"knitr",
"Gviz",
"GenomicRanges",
"survival",
"scales",
"RColorBrewer"

)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) {
  
  source("http://bioconductor.org/biocLite.R")
  biocLite(new.packages)
  
}


