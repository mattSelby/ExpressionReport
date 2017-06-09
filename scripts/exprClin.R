


exprClinicalPlot <-  function(pos,gene.of.interest, hgnc.id, clin.feat, annot, vsd, include.nos = FALSE, standard.scale = FALSE, subgroup.include = "all", clin.stats = clin.stats, output.dirname = output.dirname) {
    
  ### Get the clinical feature that is being analysed and match it to the gene.of.interest
  #
  # function to create expression vs clinical features (VSD)
  # this function is looped over so the pos is generated in the Rmd for each clinical feature to assess
  #
  # pos = which clinical feature will be used, gives a position (numeric)
  # gene.of.interest = gene of interest ENSEMBL
  # hgnc.id = the gene id created earlier
  # clin.feat = a list of the clinical features which will be compared, generated earlier
  # annot= expression annotation file
  # vsd = expression vsd file
  # include.nos = include NOS samples (TRUE)
  # standard.scale = use a standard scale (TRUE)
  # subgroup.include =  which subgroups to include but used to tell which data is INCLUDED in the plot not for comparison, all, WNT, SHH, Grp3, Grp4
  #^^^^^ case sensitive! ^^^^^
  # clin.stats = is a list of statistical comparions to be performed for each clinical feature, this can be T test or ANOVA, specified before hand in a list (please see other file)
  # output.dirname = location for the log file
  #
  
  # find the position of the clinical feature to be assessed
  
    which(pos == clin.feat) -> p
    
  # if it is the first one in the clinical features list it will create a log file and start sinking
    if (p == 1) {
      LogSinker(
        data.type = "ClinFeats", outdir = paste0(output.dirname,"/ClinicalFeats/"), start = TRUE, cont = FALSE
      )
      
    } else {
  # if it is not the first in the list it will append to the previous list
      LogSinker(
        data.type = "ClinFeats", outdir = paste0(output.dirname,"/ClinicalFeats/"), start = TRUE, cont = TRUE
      )
      
    }
    
  
  #
  # this bit evaluates the parsed text.....
  # i.e. paste0(deparse(substitute(annot)), "$", sex)
  # gives annot$sex
  # the eval(parse(text = "annot$sex")) then evaluates this to load in the data
  #
    eval(parse(text = (paste0(
      deparse(substitute(annot)), "$", clin.feat[p]
    )))) -> temp.clin
    
  # annot$Sex -> temp.clin  
    
  # which clin feat for the log  
    
    cat("\n", clin.feat[p], "\n\n")
    
  # find the statistics to be performed from the list  
    
    temp.clin.stats <-
      clin.stats[grep(clin.feat[p], names(clin.stats))]
    
    
    # sort out the names and find the match
    gsub('\\.[0-9]+$', '', rownames(vsd)) -> gene.names
    which(gene.of.interest == gene.names) -> gene.match
    
    # stop the function if no data is found
    if (length(gene.match) < 1) {
      stop("no gene match")
    }
    else if (length(gene.match) > 2) {
      stop("more than one gene matches")
    }
    
    
    ### work out which subgroups are in the analysis and if NOS is in
    
    # pull out subgroup and NMB names
    annot$Group -> subgroup
    annot$Index -> nmb.name
    
    #
    # which subgroup data to include in the analysis and whether NOS are included, also creates the titles
    #
    if (paste0(subgroup.include, collapse = "|") == "all") {
      subgroup.data -> subgroup.include
      if (include.nos == TRUE) {
        subgroup.include <- c(subgroup.include, "NOS")
      }
      
      if (is.character(hgnc.id[1,]) == "TRUE") {
        title = paste0(gene.of.interest," (",hgnc.id,")\n Expression vs ",gsub("_", " ", clin.feat[p]))
      }
      else if (is.character(hgnc.id[1,]) == "FALSE") {
        title = paste0(gene.of.interest,"\n Expression vs ", gsub("_", " ", clin.feat[p]))
      }
      
    } else {
      if (is.character(hgnc.id[1,]) == "TRUE") {
        title = paste0(
          gene.of.interest," (",hgnc.id,") expression (",paste(subgroup.include, collapse  = ", "),") vs ", 
          gsub("_", " ", clin.feat[p])
        )
      }
      else if (is.character(hgnc.id[1,]) == "FALSE") {
        title = paste0(gene.of.interest," expression (",paste(subgroup.include, collapse  =
                                                   ", "),") vs Others")
      }
      
      if (include.nos == TRUE) {
        subgroup.include <- c(subgroup.include, "NOS")
      }
      
      
    }
    
    
    # find the data to be included
    
    which(subgroup %in% subgroup.include) -> keep.index
    
    temp.clin <-  temp.clin[keep.index]
    vsd.keep <- vsd[,keep.index]
    subgroup <- subgroup[keep.index]
    nmb.name <- nmb.name[keep.index]
    ### Sort out the scale
    
    # make it a number for plotting
    as.numeric(vsd.keep[gene.match,]) -> gene.exp
    
    # standard scale?
    
    if (standard.scale == TRUE) {
      scale <- c(-1,21)
    }else{
      scale <- c(min(gene.exp - 0.5),max(gene.exp) + 0.5)
    }
    
    
    ### Remove NA/sort factor out
    ### relapse is dealt with differently, add other in
    if (clin.feat[p] == "Relapse") {
      temp.clin <- as.factor(gsub("^$","Other",temp.clin))
      temp.clin <-
        factor(temp.clin,levels = c("Relapse","Progression", "Other"))
      table(temp.clin) -> numbers
      
    } else {
      ### Sort out the rest if not relapse
      # get rid of any whitespace with NA
      # then strip NA values out
      
      temp.clin <- as.factor(gsub("^$","NA",temp.clin))
      which(temp.clin == "NA") -> na.index
      num.na <- length(na.index)
      temp.clin <- temp.clin[-na.index]
      gene.exp <- gene.exp[-na.index]
      temp.clin <- factor(temp.clin)
      table(temp.clin) -> numbers
      subgroup <- subgroup[-na.index]
      nmb.name <- nmb.name[-na.index]
      
    }
    
    
    #### axis size
    
    # use the length of clinical feature to get axis size
    
    if (length(levels(temp.clin)) > 5) {
      ax.size = .6
    } else {
      ax.size = 1
    }
    # function to make a lovely salmon to firebrick gradient...can be changed
    colfunc <- colorRampPalette(clinical.colours)
    
    ### title and colour
    
    # work out how many of each sample....
    names <- paste0(levels(temp.clin), " (n=", numbers, ")")
    
    # make colours
    cols <- colfunc(length(levels(temp.clin)))
    
    # numbers for log
    cat("\nNumbers: ", names, "\n\n")
    
    ## plotting and median line/text
    par(mar = c(8, 6, 5, 6), las = 1)
    par(xpd = FALSE)
    boxplot(
      gene.exp ~ temp.clin, col = cols, ylim = scale, ylab = "Expression (VSD)", names = names, main = title, xaxt = 'n'
    )
    axis(
      side = 1,at = c(1:length(levels(temp.clin))),labels = names, cex.axis =  ax.size
    )
    abline(h = median(gene.exp), lty = 2)
    # write in the margins
    par(xpd = TRUE)
    
    # median text written
    text(par("usr")[2],median(gene.exp),"Median", cex = 1, adj = -.1)
    
    #
    # this loop is for working out the stats the user put in
    # the k loop goes round as many times as there are seperate test to perform
    # e.g male, female, other
    # tests to be compared:
    # male, female vs other
    # male vs female
    # male vs female vs other
    # would loop over 3 times
    #
    
    # temp.clin.stats <- list(sex = c("M", "F"))
    # k=1
    # j=11
    
    for (k in 1:length(temp.clin.stats)) {
      # create output lists
      temp.clin.forStats <- rep(NA, length(temp.clin))
      names.forStats <- list()
      
    #
    # this little loop replaces the different selected options from the clinical features with numbers, 
    # e.g. M= 1 F=2
    # and also creates the names used in the stats testing later  
    # it leaves and NAs out  
      
      for (j in 1:length(temp.clin.stats[[k]])) {
        temp.clin.forStats[grep(paste0(temp.clin.stats[[k]][[j]], collapse = "|"), temp.clin)] <-   j
        names.forStats[[j]] <- paste0(temp.clin.stats[[k]][[j]], collapse = "|")
        
      }
      
      # remove any of the input data that is an NA, e.g. any data not to be include would have an NA
      gene.exp.forStats <- gene.exp[!is.na(temp.clin.forStats)]
      temp.clin.forStats <- temp.clin.forStats[!is.na(temp.clin.forStats)]
      
      # work out the numbers and for each componenet
      
      numbers.forStats <- table(temp.clin.forStats)
      names.forStats <- paste0(names.forStats, " (n=", numbers.forStats, ")")
      
      # write that to the log file
      cat("\nNumbers to test: ", paste0(names.forStats, collapse = " vs. "), "\n\n")
      
      ### t test for only two variables
      # if there are only two variables:
      
      if (length(temp.clin.stats[[k]]) == 2) {
        ### make sure there are sufficient numbers for analysis
        
        # this makes sure there are enough to do stat testing
        if (!any(numbers.forStats <= 1)) {
          
          try(t.test(gene.exp.forStats ~ temp.clin.forStats, silent=TRUE) -> stat)
          
          
          if (is(stat, "try-error")) {stat <-NA }
          
          
          # log the stats test
          cat("\n T test performed \n\n")
          
          
          if(!all(is.na(stat))){
          (print(stat))
          
          # pull out the values
          stat[[2]] -> df
          stat[[3]] -> p.val
          
          ### work out text and plotting of that text
          
          if (p.val < 0.001 & !is.na(p.val)) {
            p.text <- "p < 0.001"
          } else {
            p.text <- paste("p =" , round(p.val,3))
          }
          
          f.text <- paste("DF =",signif(df,3))
          
          # write on the chart
          mtext(
            paste0(
              "T Test for ", paste0(names.forStats, collapse = " vs. "),  "; ", f.text,", " ,p.text
            ),
            at = par("usr")[1],
            line = -33 + 1 - k, 
            adj = 0, cex = .8
          )
          } 
        } else {
          
          # this tells you there aren't enough to test in the margin 
          mtext(
            paste0(
              "T Test for ", paste0(names.forStats, collapse = " vs. "),  "; ", "Too few to test"
            ), at = par("usr")[1], line = -33 + 1 - k, adj = 0, cex = .8
          )
          
          
          
        }
        
      # this is for multiple samples  
        
      }  else if (length(temp.clin.stats[[k]]) > 2) {
        ### sames as above but for multiple variables, ANOVA
        
        # again check enough numbers are present
        if (!any(numbers.forStats <= 1)) {
          # do the ANOVA
          try(aov(gene.exp.forStats ~ temp.clin.forStats) -> stat)
          
          if (is(stat, "try-error")) {stat <-NA }
          if(!all(is.na(stat))){
          try(summary(stat) -> stat.summary)
          
          # output to log file
          cat("\n ANOVA performed \n\n")
          
          try(print(stat.summary))
          
          # pull out the stats and write on  the chart
          stat.summary[[1]]$"F value"[1] -> F
          stat.summary[[1]]$"Pr(>F)"[1] -> p.val
          
          
          ### work out text and plotting of that text
          
          if (p.val < 0.001 & !is.na(p.val)) {
            p.text <- "p < 0.001"
          } else {
            p.text <- paste("p =" , round(p.val,3))
          }
          
          
          f.text <- paste("F =",round(F,3))
          
          mtext(
            paste0(
              "ANOVA for ", paste0(names.forStats, collapse = " vs. "),  "; ", f.text,", " ,p.text
            ), at = par("usr")[1], line = -33 + 1 - k, adj = 0, cex = .8
          )
          }
        } else {
          # tell the user there are too few to test
          mtext(
            paste0(
              "ANOVA for ", paste0(names.forStats, collapse = " vs. "),  "; ", "Too few to test"
            ), at = par("usr")[1], line = -33 + 1 - k, adj = 0, cex = .8
          )
          
        }
        
      }
      
      
    }
    
    # this is to make an output csv which gets saved later after output
    temp.clin.output <- rbind(as.character(temp.clin), as.character(subgroup), gene.exp)
    rownames(temp.clin.output) <- c(clin.feat[[p]], "Subgroup", "Expression")
    colnames(temp.clin.output) <- nmb.name
    # turn the sink off and export
    LogSinker(start = FALSE)
    return(temp.clin.output)
    
  }
