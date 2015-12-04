

LogSinker <-  function(data.type, outdir = outdir, start = TRUE, cont = FALSE) {
  # script to create the log files that are output into the out dirs
  # data.type = type of data... e.g. Expression
  # outdir = the output dir for the log file, this is often hardcoded in the scripts for ease as no one should be changing them
  # start = is this the start (TRUE) of the log file?
  # cont = do you want to append to a preexisting log file (TRUE)?
  #
    
    if (start == TRUE) {
      if (cont == FALSE) {
        #start the sink by creating the log
        sink(paste0(outdir,"LogFile", data.type, ".Rout"))
        
        # print the R version
        print(as.matrix(R.Version()))
        
        # with spaces print the below details
        cat("\n")
        cat("Gene of Interest = ", gene.of.interest, "\n")
        cat("HGNC ID = ", hgnc.id[[1]], "\n")
        cat("Group of Interest = ", group.of.interest, "\n")
        cat("Included Subgroups = ", subgroup.include, "\n")
        cat("Standard Scale = ", standard.scale, "\n")
        cat("Include NOS = ", include.nos, "\n")
        cat("\n")
        
      } else {
        # or add it to a preexisting one
        sink(paste0(outdir,"LogFile", data.type, ".Rout"), append = TRUE)
        
      }
      
    } else {
      # turn the sink off
      sink()
      
    }
  }
