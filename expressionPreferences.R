#
#R version 3.4.0 (2017-04-21)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 16.04.2 LTS
#
#Matrix products: default
#BLAS: /usr/lib/libblas/libblas.so.3.6.0
#LAPACK: /usr/lib/lapack/liblapack.so.3.6.0
#
#locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
#[5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       #
#
#attached base packages:
#  [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
 # [1] RColorBrewer_1.1-2         scales_0.4.1               survival_2.41-3            Gviz_1.20.0               
#[5] rtracklayer_1.36.0         minfi_1.22.1               bumphunter_1.16.0          locfit_1.5-9.1            
#[9] iterators_1.0.8            foreach_1.4.3              Biostrings_2.44.0          XVector_0.16.0            
#[13] SummarizedExperiment_1.6.1 DelayedArray_0.2.2         matrixStats_0.52.2         GenomicRanges_1.28.1      
#[17] GenomeInfoDb_1.12.0        reshape2_1.4.2             biomaRt_2.32.0             org.Hs.eg.db_3.4.1        
#[21] AnnotationDbi_1.38.0       IRanges_2.10.0             S4Vectors_0.14.0           Biobase_2.36.2            
#[25] BiocGenerics_0.22.0        markdown_0.8               rmarkdown_1.5              xtable_1.8-2              
#[29] knitr_1.15.1              
#
#loaded via a namespace (and not attached):
#  [1] nlme_3.1-131                  ProtGenerics_1.8.0            bitops_1.0-6                 
#[4] httr_1.2.1                    rprojroot_1.2                 tools_3.4.0                  
#[7] backports_1.0.5               doRNG_1.6.6                   nor1mix_1.2-2                
#[10] R6_2.2.0                      rpart_4.1-11                  Hmisc_4.0-3                  
#[13] DBI_0.6-1                     lazyeval_0.2.0                colorspace_1.3-2             
#[16] nnet_7.3-12                   gridExtra_2.2.1               base64_2.0                   
#[19] compiler_3.4.0                preprocessCore_1.38.1         htmlTable_1.9                
#[22] pkgmaker_0.22                 checkmate_1.8.2               genefilter_1.58.1            
#[25] quadprog_1.5-5                stringr_1.2.0                 digest_0.6.12                
#[28] Rsamtools_1.28.0              foreign_0.8-67                illuminaio_0.18.0            
#[31] siggenes_1.50.0               GEOquery_2.42.0               dichromat_2.0-0              
#[34] base64enc_0.1-3               htmltools_0.3.5               ensembldb_2.0.1              
#[37] BSgenome_1.44.0               limma_3.32.2                  htmlwidgets_0.8              
#[40] RSQLite_1.1-2                 BiocInstaller_1.26.0          shiny_1.0.0                  
#[43] mclust_5.2.3                  BiocParallel_1.10.1           acepack_1.4.1                
#[46] VariantAnnotation_1.22.0      RCurl_1.95-4.8                magrittr_1.5                 
#[49] GenomeInfoDbData_0.99.0       Formula_1.2-1                 Matrix_1.2-10                
#[52] Rcpp_0.12.10                  munsell_0.4.3                 stringi_1.1.5                
#[55] yaml_2.1.14                   MASS_7.3-47                   zlibbioc_1.22.0              
#[58] AnnotationHub_2.8.1           plyr_1.8.4                    lattice_0.20-35              
#[61] splines_3.4.0                 multtest_2.32.0               GenomicFeatures_1.28.0       
#[64] annotate_1.54.0               beanplot_1.2                  rngtools_1.2.4               
#[67] codetools_0.2-15              XML_3.98-1.7                  evaluate_0.10                
#[70] biovizBase_1.24.0             latticeExtra_0.6-28           data.table_1.10.4            
#[73] httpuv_1.3.3                  gtable_0.2.0                  openssl_0.9.6                
#[76] reshape_0.8.6                 ggplot2_2.2.1                 mime_0.5                     
#[79] AnnotationFilter_1.0.0        tibble_1.3.0                  GenomicAlignments_1.12.0     
#[82] registry_0.3                  memoise_1.1.0                 cluster_2.0.6                
#[85] interactiveDisplayBase_1.14.0









#
# The preferences file is used to specify input/plotting data to the report
# script, it is set to a default that produces reports for all subgroups but can
# be altered for a more bespoke report
#
# There are three sections below:
#

#
# DATA TO INCLUDE IN THE REPORT
#
{
# These preferences allow the user to specify which data to compare/include
# within the report
#


#
# Subgroup Information
#
# In the case of medulloblastoma this is set for the four subgroups but these
# can be changed to suit the data that has been input
# 


#
# Subgroups sets the list of subgroups, these must be written exactly the same
# as in the annotation file, the order of these subgroups also determines the
# order you will plot them in on the graphs
#

subgroup.data <- c("WNT", "SHH", "Grp3", "Grp4")




#
# Subgroup Include is used to remove specific data from the analysis
#
# By default all subgroups ("all") are included in the analysis but data can be removed
#
# e.g. "WNT" would only show WNT data in the report
# e.g. c("SHH", Grp3", "Grp4") would remove the WNT data from the report
#
#  options : c("all", "WNT", "SHH", "Grp3", "Grp4")


subgroup.include <- "all"

#
# group of interest: this variable is used to specify whether there is a group of
# interest
# For example: if you wanted to look across all subgroups (which is the
# default) group.of.interest = "all"
# An ANOVA would be performed
#
# Group of Interest can be used to compare individual subgroups vs. the rest of the data set
# "Grp3" would set all other subgroups to Other and compare that to Grp3
#
# You can also compare multiple subgroups e.g c("Grp3", "Grp4") would compare
# Grp3/4 vs other subgroups (in this case "SHH" and "WNT")
#
# options : c("all", "WNT", "SHH", "Grp3", "Grp4")
#

group.of.interest <- "all"

#
# In the medullo set there are NOS (No Official Subgroup) samples, these by
# default are removed  but can be left in the analysis if required
#
# options = TRUE/FALSE
#

include.nos <- FALSE

#
# surv.plots gives a list of the individual survival plots/ COX models to plot
# default is to plot all data (c("all", "WNT", "SHH", Grp3", "Grp4") together then individual plots for each subgroup
#
# can choose to use one subgroup...e.g. list("Grp3") will only give the Grp3 plot
# or can group together
# e.g. list("all", c("Grp3", "Grp4"))
# would give the plot with all data on then a second plot with just Grp3 and
# Grp4
#
# options : c("all", "WNT", "SHH", "Grp3", "Grp4")
#

surv.plots <- list("all", "WNT", "SHH", "Grp3", "Grp4")

#
# Clin Feat specifies which of the clinical features from the annotation files
# should be plotted
#
# There is a default list provided below, which also corresponds to the clinical
# stats list
#
# The helper function below can be used to look at all the fields from the annotation
# file and then select which ones are carried forward into the report.
#
# To run this script remove the hashes, first run the "script.location" line of
# code and the "db" line of code in the DATA INPUT section below
#
# Following this run the source code and then the clinFeatHelperFunction, this
# will give a list of all the columns in the annotation from which the user can
# choose from
#
# If a change from the default is required use the exact same format from the
# list and change the clin.feat below
#
# CLIN FEAT HELPER CODE
#
# source(paste0(script.location, "clinFeatHelperFunction.R"))
# clinFeatHelperFunction(db)

clin.feat <- c(
  "Sex",                  
  "Age",
  "Stage",
  "Resection", 
  "Status_follow_up",
  "Relapse",           
  "M_status_at_relapse", 
  "Beta_catenin_IHC_result", 
  "Beta_catenin_Sequencing_result",
  "TP53_IHC_result",
  "TP53_Sequencing_result",     
  "MYC_FISH_result",             
  "MYCN_FISH_result",     
  "TERT_Sequencing_result", 
  "TERT_Taqman_result",
  "Final_pathology"
)    


#
# CLINICAL FEATURE NAMES SHOULD BE USED (EXACTLY!!)
#
# The clinical stats list is used to perform statistical analyses on the
# clinical plots which are produced, there are a default set of plots (as below)
# but with some alterations this list can be used to interrogate the data in
# many different ways
#
# For one analysis per feature the name of the Clinical Feature (exactly) is used in the format below:
#
# Sex = list("F", "M")
#
# T test between female and male
#
# OR
#
# Stage = list(c("M0",   "M0/1", "M0/2", "M0/3", "M1"), c("M2",  "M3", "M4")),
#
# The groupings within the list are also have a t test but between all in those groups
# e.g "M0",   "M0/1", "M0/2", "M0/3", "M1" vs. "M2",  "M3", "M4"
#
# If multiple analysis of the same feature are required the following format is used:
#
# Stage_1 = list(c("M0",   "M0/1", "M0/2", "M0/3", "M1"), c("M2",  "M3", "M4"))
#
# Stage_1 would do a t test between the two groups
# e.g "M0",   "M0/1", "M0/2", "M0/3", "M1" vs. "M2",  "M3", "M4"
#
# Stage_2 = list("M0",   "M0/1", "M0/2", "M0/3", "M1", "M2",  "M3", "M4")
#
# Stage_2 would do an ANOVA across all entries
# e.g "M0" vs.  "M0/1" vs. "M0/2" vs. "M0/3" vs. "M1" vs. "M2" vs.  "M3" vs. "M4"
#
# Stage_3 = list(c("M0",   "M0/1", "M0/2", "M0/3"), c("M1", "M2"),  c("M3", "M4"))
#
# Stage_3 would do an ANOVA across the three groups
# e.g "M0",   "M0/1", "M0/2", "M0/3" vs. "M1" vs. "M2",  "M3", "M4"
#


# The names below each list entry are all the options for that particular
# clinical feature

clin.stats <- list(
  #
  Sex = list("F","M"),
  #
  #  "F", "M"
  #
  Age = list("adult", "child", "infant"),
  #
  #  "adult", "child", "infant"
  #
  Stage_1 = list(c("M0",   "M0/1", "M0/2", "M0/3", "M1"), c("M2",  "M3", "M4")),
  Stage_2 = list(c("M0",   "M0/1", "M0/2", "M0/3"), c("M1","M2",  "M3", "M4")),
  #
  #  "M0",  "M0/1", "M1",  "M2",  "M2/3", "M3", "M4"
  #
  Resection = list("GTR", "STR"),
  #
  #  "Biopsy", "CR", "GTR", "STR"
  #
  Status_follow_up = list(c("AWD", "DOD"), c("ADF","DOOC")),
  #
  #  "ADF", "AWD", "DOD", "DOOC"
  #
  Relapse_1 = list("Other", "Relapse"),
  Relapse_2 = list("Progression", "Relapse"),
  #
  #   If others is selected then takes into account entires with nothing in AND progression
  #
  #  "Progression", "Relapse", "Other"
  #
  M_status_at_relapse = list(c("M0",  "M0/1", "M0/2", "M0/3", "M1"), c("M2",  "M3", "M4")),
  #
  #  "M0",  "M0/1", "M1",  "M2",  "M2/3", "M3", "M4"
  #
  Beta_catenin_IHC_result  = list("Negative", "Positive"),
  #
  #  "Negative", "Positive"
  #
  Beta_catenin_Sequencing_result = list("MUT", "WT"),
  #
  #  "MUT", "WT"
  #
  TP53_IHC_result = list("Negative", "Positive"),
  #
  #  "Negative", "Positive"
  #
  TP53_Sequencing_result = list("MUT", "WT"),
  #
  #  ""MUT", "WT"
  #
  MYC_FISH_result = list("Amplified", "Balanced"),
  #
  #  "Amplified", "Balanced"
  #
  MYCN_FISH_result = list("Amplified", "Balanced"),
  #
  #  "Amplified", "Balanced"
  #
  TERT_Sequencing_result = list("Mut", "WT"),
  #
  #  "MUT", "WT"
  #
  TERT_Taqman_result = list("Mut", "WT"),
  #
  #  "MUT", "WT"
  #
  Final_pathology_1 = list("CLA", "LCA"),
  Final_pathology_2 = list("CLA", "DN")
  #
  #  "CLA", "DN", "LCA", "MBEN", "NOS", "Other"
  #
)

## The bin allows for customisability of the KM, options are:
# "median", this is default
# "mean"
# "quart.100", 100% quartile
# "quart.75", 75% quartile
# "quart.50", 50% quartile
# "quart.25", 25% quartile
# "quart.0", 0% quartile
# 
# any number, this takes the data as a percentage
#

#bin <- c("median", "mean","quart.75", "quart.50", "quart.25")
bin = "median"
}

#
# Not sure about adding in the ability to split the data into more than one subset
# or to remove some data e.g. <25% and >75%...will give it a go
#
# going to try, not sure how to deifne if it is a quartile/tertile etc. will start with just quartile
# split indicates that the data is to be split into 3 sets
# 
# bin = "split_25_75"
# outlier indicates that the data will be below and above the two values
#
# bin = "outlier_25_75"



#
# SCALES, COLOURS AND LIMITS
#
{
# This section allows the user to change the plotting scales?! and also limits
# for a couple of scripts
#


#
# Cutoff is used in the fusion script, this number specifies how many
# spanning/junction reads are required to be classed as a fusion, default is set
# to 2
#
cutoff <- 2

#
# Flag Limit is used in the survival analysis, if there is more than the flag
# limit (default = 75%) of data missing an error message is printed
#

flag.limit <- .75

# EXPERIMENTAL FEATURE
#
# Do you want to use the same scale across all the plots, doesn't work as there
# are both VSD and FPKM which can vary
#
# Will probably remove this
#

standard.scale <- FALSE


#
# Colour Information
#
# All colours have defaults in the report but if changes are required the lists
# below can be altered
#

#
# Subgroups colours is used to determine plotting colours later in the scripts,
# these are default to a lovely pastel palette but can be changed here. These
# also contain colours for Other and NOS which are currently set to grey
#
# This is a named list, the subgroup must always be on the right, then an
# equals, then the colour for that subgroup

subgroup.colours <- list(
  
  "WNT" = "steelblue2",
  "SHH" = "tomato3",
  "Grp3" = "gold1",
  "Grp4" = "darkolivegreen3",
  "Other" = "grey",
  "NOS" = "grey",
  "CB" = "grey"
  
)

# 
# Clinical colours are used to set the two colours which generate a colour scale
# for each of the clincal feature plots
#

clinical.colours <- c("light salmon","firebrick")

#
# KM colours (kaplan meier) set the colours for the two lines on the KM curves
#

km.colours <- c("red", "blue", "green")

}


#
# DATA INPUT
#
{
# This section holds locations of data files to be imported, this is used to
# update data sets and gives folder locations for data to be imported.
# It also contains the location of the folder containing plotting scripts


# This folder contains all of the plotting scripts required to generate the
# expression report, currently this will need changing when you first import the
# report scripts
script.location <- "/home/shelby/ExpressionReport/scripts/"

# This folder has the svg fusion output from SOAP to extract if there are
# fusions found
#
# Currently hardcoded to Dan's folder but this needs centralising then will be
# redundant

SOAPdir = "/home/dan/SOAP/output"


# This folder contains all of the various subfolders with the data to plot in
# the report, to update the data a second folder will be generated with the data
# and data in the same format. This will need changing as and when updates are
# rolled out
#
# NOTE
#
# IN THE LOAD DATA SCRIPT THERE ARE SOME FILES THAT NEED GERENATING BUT ONLY
# ONCE, WHEN AN UPDATE IS APPLIED THE LAOD DATA SCIRPT WILL NEED ALTERING TO
# GENERATE THESE DATA STRUCTURES
#

db <- "/home/data/pbt/RNASeq/report_data/MB_20_April_2015"

}

convertIDs<- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg(ifMultiple)
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}
