aseTableCreator <- function(gene.of.interest, output.dirname = output.dirname, ase.data = ase.data, subgroup.include = "all", include.nos = FALSE) {
    
  # the data provided for the ase has pre processsed and non processed ase data
  # I chose to use the non processed and do it myself
  #
  # Function creates the table of relative allele specific expression for the gene of interest
  
  # gene.of.interest = gene of interest (ENSEMBL)
  # output.dirname = location for the log file
  # ase.data = the data provided for the allele specific expression
  # subgroup.include =  which subgroups to include but used to tell which data is INCLUDED in the plot not for comparison, all, WNT, SHH, Grp3, Grp4
  #^^^^^ case sensitive! ^^^^^
  # include.nos = include NOS samples (TRUE)
  #
  #
  
  
    #[[1]]    res.mb          res.mb is results matrix TRUE/FALSE/NA
    #[[2]]    subgroup        subgroup is vector of subgroups
    
    # gene.of.interestNG TO CALCULATE MY OWN FROM HERE
    
    #[[3]]    num.ase.mb      num.ase.mb is a vector of the number of TRUE cases per gene
    #[[4]]    num.test.mb     num.test.mb is the number of TRUE + FALSE cases per gene
    #[[5]]    perc.ase.mb     perc.ase.mb is a vector of the percentage of TRUE cases over the total number of determined cases
    
    
    # turn on the sink
    
    LogSinker(
      data.type = "ASE", outdir = paste0(output.dirname,"/Expression/ASE/"), start = TRUE, cont = FALSE
    )
    
    
    
    # do we keep nos?
    if (include.nos == FALSE) {
      # pull out the subgroup and remove NOS
      subgroup <- ase.data[[2]]
      nos.index <- which(subgroup == "NOS")
      subgroup <- subgroup[-nos.index]
      
      # get the ase for the gene of interest
      gene.of.interest.ase <- ase.data[[1]][which(rownames(ase.data[[1]]) == gene.of.interest),-nos.index]
      
      # find all the True/All and Na cases from the ASE
      TrueCases <- sum(gene.of.interest.ase, na.rm = TRUE)
      AllCases <-  sum(!is.na(gene.of.interest.ase))
      NaCases <- sum(is.na(gene.of.interest.ase))
      
      # work our the perecentage of true cases
      PercCases <- TrueCases / AllCases
      
      # trim down the subgroup and ase for the gene of interest
      subgroup <- subgroup[!is.na(gene.of.interest.ase)]
      gene.of.interest.ase <- gene.of.interest.ase[!is.na(gene.of.interest.ase)]
      
      
      # work out which subgroup to include in the analysis
      if (subgroup.include == "all") {
        c("WNT","SHH","Grp3","Grp4") -> temp.subgroup.include
        
      } else {
        subgroup.include -> temp.subgroup.include
        
      }
    } else {
      
      # as above but not removing the NOS samples
      
      subgroup <- ase.data[[2]]
      
      gene.of.interest.ase <- ase.data[[1]][which(rownames(ase.data[[1]]) == gene.of.interest),]
      
      TrueCases <- sum(gene.of.interest.ase, na.rm = TRUE)
      AllCases <-  sum(!is.na(gene.of.interest.ase))
      NaCases <- sum(is.na(gene.of.interest.ase))
      PercCases <- TrueCases / AllCases
      
      subgroup <- subgroup[!is.na(gene.of.interest.ase)]
      gene.of.interest.ase <- gene.of.interest.ase[!is.na(gene.of.interest.ase)]
      
      if (subgroup.include == "all") {
        c("WNT","SHH","Grp3","Grp4","NOS") -> temp.subgroup.include
      } else {
        subgroup.include -> temp.subgroup.include
        
      }
      
    }
    
  # turn this into a factor to keep the order right in the plot
    subgroup <- factor(subgroup, levels = temp.subgroup.include)
    # subgroup keep index
    keep.index <- which(subgroup %in% temp.subgroup.include)
    
    # work out the relative proportions for the subgroups
    sg.prop <- as.matrix(prop.table(table(subgroup)))
    
    # get the number positions
    num.positions <- length(rownames(sg.prop))
    
    # reverse to get WNT first, then work out the ratios
    proportions <- as.matrix(rev(prop.table(table(gene.of.interest.ase))))
    
    # this set of ifs make the colour scheme
    # hard coded but can be manually changed if pastels aren't your thing
    sbgrp.col <- c()
    if (any(levels(subgroup) == "Other")) {
      sbgrp.col  <- "grey"
    }
    if (any(levels(subgroup) == "WNT")) {
      sbgrp.col  <- c(sbgrp.col,"steelblue2")
    }
    if (any(levels(subgroup) == "SHH")) {
      sbgrp.col <- c(sbgrp.col,"tomato3")
    }
    if (any(levels(subgroup) == "Grp3")) {
      sbgrp.col <- c(sbgrp.col,"gold1")
    }
    if (any(levels(subgroup) == "Grp4")) {
      sbgrp.col <- c(sbgrp.col,"darkolivegreen1")
    }
    if (any(levels(subgroup) == "NOS")) {
      sbgrp.col <- c(sbgrp.col,"grey")
    }
    
    # create the two output vectors
    ratios <- c()
    cols2 <- c()
    # move over the number of subgroups in sequence gene.of.interestng two at a time
    for (i in seq(1,length(levels(subgroup)) * 2,2)) {
      
      #  looks a bit complicated but take the first bit and work out the ratio of e.g. subgroup
      ratios[i] <- (prop.table(table(subgroup)) * tapply(gene.of.interest.ase, subgroup, mean))[(i + 1) /2]
      # give it a colour
      cols2[i] <- sbgrp.col[(i + 1) / 2]
      # make the other ratio to plot
      ratios[i + 1] <- (prop.table(table(subgroup)) * (1 - tapply(gene.of.interest.ase, subgroup, mean)))[(i + 1) / 2]
      # give it an alphaed colour
      cols2[i + 1] <- adjustcolor(col = sbgrp.col[(i + 1) / 2], alpha.f = 0.2)
      
    }
    
    # turn NA into 0
    ratios[is.na(ratios)] <- 0
    
    if(sum(ratios) >0){
    # get plotting
    par(mar = c(0,8,12,10))
    
    # make the overall plot of TRUE/FALSE
    barplot(
      proportions, las = 1,
      axes = FALSE, horiz = TRUE, width = .5, ylab = "", col = c("#606060","grey"), ylim = c(0,5), space =
        9
    )
    
    # put on a smaller subgroup bar
    barplot(
      sg.prop, add = T, axes = F, horiz = TRUE, width = 0.25, ylab = "", col = sbgrp.col, space = 15
    )
    
    # put on the ratios for each of the subgroup of TRUE/FALSE
    barplot(
      as.matrix(ratios), add = T, axes = F, horiz = TRUE, width = 0.5, ylab = "", col = cols2, space = 6.5
    )
    
    # plot an axis...
    axis(
      3, at = c(0,0.2,0.4,0.6,0.8,1),labels = c("0%","20%","40%","60%","80%","100%")
    )
    
    # give numbers of true cases etc.
    mtext(
      text = paste0(
        "True Cases: ", TrueCases, " ; False Cases: ", AllCases - TrueCases," ; NA: ", NaCases
      )
      , at = c(0,0), line = -6.5, cex = 0.8, adj = 0
    )
    
    # more margin stuff, this is hard coded for knitr so would look bad out of the .Rmd
    mtext(
      text = "Subgroup", side = 2, las = 2, adj = 1, padj = -19, cex = 0.8
    )
    mtext(
      text = "True/False", side = 2, las = 2, adj = 1, padj = -13.5, cex = 0.8
    )
    # get the lgened bits ready
    leg.lab <- c("TRUE", "FALSE", "\n", temp.subgroup.include)
    leg.col <- c("#606060", "grey")
    leg.ang <- c(45, 135)
    
    # plot the legends
    legend(
      x = 1.1, y = par("usr")[4] - .3, legend = c("TRUE", "FALSE"),
      fill =  leg.col, xpd = TRUE, cex = .8,
      title = "Proportion of \nTrue/False Calls\n from ASE", bty = "n", xjust = .5
    )
    
    
    legend(
      x = 1.1, y = par("usr")[4] - 1.25, legend = temp.subgroup.include,
      fill =  sbgrp.col, xpd = TRUE,  cex = .8,
      title = "Subgroup", bty = "n", xjust = .5
    )
    
    }
    # get a count of the various ratios for output
    # this is transposed for aesthetic
    sbgrp.counts <-
      as.character(t(as.matrix((unlist(
        lapply(X = split(factor(
          gene.of.interest.ase, levels = c("TRUE", "FALSE")
        ), subgroup), FUN = table)
      )))))
    
    # turn the ratios into % and table up
    sbgrp.percs <-
      unlist(lapply(
        X = split(factor(gene.of.interest.ase, levels = c("TRUE", "FALSE")), subgroup), FUN = function(x) {
          prop.table(table(x))
        }
      ))
    sbgrp.percs[is.na(sbgrp.percs)] <- 0
    sbgrp.counts <- rbind(sbgrp.counts, round((sbgrp.percs * 100),2))
    
    colnames(sbgrp.counts) <- gsub("\\.", " ", colnames(sbgrp.counts))
    rownames(sbgrp.counts) <- c("Count", "Percentage (%)")
    
    #sbgrp.counts[1,] <- round(sbgrp.counts[1,])
    
    # put it all in the output file
    ASE.output <- sbgrp.counts
    
    # this is for the log file
    cat("\nTrue Cases:\n")
    cat(TrueCases)
    cat("\nAll Cases:\n")
    cat(AllCases)
    cat("\nNA Cases:\n")
    cat(NaCases)
    cat("\nPercentage Cases:\n")
    cat(PercCases)
    cat("\nCount Table:\n\n")
    cat(sbgrp.counts)
    
    # stop the sink and return the other bits
    LogSinker(start = FALSE)
    
    return(ASE.output)
    
    
  }
