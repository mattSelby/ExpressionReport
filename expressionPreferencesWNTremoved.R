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


subgroup.include <-  c("SHH", "Grp3", "Grp4")

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

surv.plots <- list("all", "SHH", "Grp3", "Grp4")

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

bin <- "median"

}



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

km.colours <- c("red", "blue")

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
