# list of packages required to run the expression report

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

# check against installed packages to see if any need to be installed before running

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) {
  
  source("http://bioconductor.org/biocLite.R")
  biocLite(new.packages)
  
} else {
  
  message("\nAll required packages installed\n")
}

rm(list.of.packages)
rm(new.packages)





