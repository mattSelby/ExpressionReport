exprSubgrpPlot <-  function(gene.of.interest, group.of.interest = "all", hgnc.id, annot, vsd, subgroup.include = "all", include.nos = FALSE, standard.scale = FALSE, output.dirname) {
    
  #
  # function to create expression vs subgroup plots (VSD)
  # gene.of.interest = gene of interest ENSEMBL
  # group.of.interest = this can either be "all" for all subgroups or one/ a combination of others (WNT, SHH, Grp3, Grp4) and is used to compare them
  # e.g. WNT vs others, Grp3 + Grp4 vs others
  # this is also controlled by subgroup include a bit, if WNT isn't included below then SHH chosen here it would be the equivilant
  # of SHH vs Grp3 and Grp 4 
  # ^^^^^ case sensitive! ^^^^^
  # hgnc.id = the gene id created earlier
  # annot= expression annotation file
  # vsd = expression vsd file
  # subgroup.include =  differs from group.of.interest but used to tell which data is INCLUDED in the plot not for comparison, same options as above
  #^^^^^ case sensitive! ^^^^^
  # include.nos = include NOS samples (TRUE)
  # standard.scale = use a standard scale (TRUE)
  # output.dirname = location for the log file
  #
    
  # start the logSinker function with start location and out dir
    LogSinker(
      data.type = "Expression", outdir = paste0(output.dirname,"/Expression/"), start = TRUE, cont = FALSE
    )
    
  # pull out subgroup and NMB names
    annot$Group -> subgroup
    annot$Index -> nmb.name
    
    # find the overlaps between NMB and the VSD columns
    which(nmb.name %in% colnames(vsd)) -> nmb.index
    which(colnames(vsd) %in% nmb.name) -> vsd.index
    
    # trim off any .1 from the end of ensembl IDs
    gsub('\\.[0-9]+$', '', rownames(vsd)) -> gene.names
    # find the gene match from the data set and get the expression for that
    which(gene.of.interest == gene.names) -> gene.match
    as.numeric(vsd[gene.match,vsd.index]) -> gene.exp
    
    # stop the function if no data is found
    if (length(gene.match) < 1) {
      stop("no gene match")
    }else if (length(gene.match) > 2) {
      stop("more than one gene matches")
    }
    
    # subset NMB/subgroup by nmb.index to remvoe non included samples
    subgroup <- subgroup[nmb.index]
    nmb.name <- nmb.name[nmb.index]
    
    
    ### sort out which subgroups are included
    
    # collapse the subgroups with a | to search
    if (!paste0(subgroup.include, collapse = "|") == "all") {
      # are NOS included?
      if (include.nos == TRUE) {
        subgroup.include <- c(subgroup.include, "NOS")
      }
      # create a keep.index for samples to be included and subset all required data by that
      which(subgroup %in% subgroup.include) -> keep.index
      
      subgroup <- subgroup[keep.index]
      subgroup <- factor(subgroup, levels = subgroup.include)
      gene.exp <-  gene.exp[keep.index]
      nmb.name <- nmb.name[keep.index]
      
      # create a list of the included subgroups if it was "all
    } else if (include.nos == TRUE) {
      subgroup.include <- c("WNT","SHH","Grp3","Grp4", "NOS")
    } else {
      subgroup.include <- c("WNT","SHH","Grp3","Grp4")
    }
    
    ### sort out which the group of interest are
    
    if (!any(group.of.interest == "all")) {
      # find which groups have been chosen to compare and compare to "other"
      groups <- c("WNT", "SHH", "Grp3", "Grp4")
      groups <-
        groups[grep(paste0(group.of.interest, collapse = "|"), groups, invert = TRUE)]
      groups <- paste0(groups, collapse = "|")
      subgroup <- gsub(groups,"Other",subgroup)
      
      # whether to include nos?
      if (include.nos == FALSE) {
        which(subgroup == "NOS") -> nos.index
        subgroup <- subgroup[-nos.index]
        gene.exp <- gene.exp[-nos.index]
        nmb.name <- nmb.name[-nos.index]
        #changeing to factor to sort out the plotting order
        subgroup <- factor(subgroup,levels = c("Other",group.of.interest))
      } else {
        subgroup <- factor(subgroup,levels = c("Other",group.of.interest, "NOS"))
      }
      
      #### sort the title out
      # create the titles here
      if (is.character(hgnc.id[1,]) == "TRUE") {
        title = paste0(gene.of.interest," (",hgnc.id,") expression (",paste(group.of.interest, collapse  =
                                                                 ", "),") vs Others")
      }
      else if (is.character(hgnc.id[1,]) == "FALSE") {
        title = paste0(gene.of.interest," expression (",paste(group.of.interest, collapse  = ", "),") vs Others")
      }
      
      
    } else {
      #### this bit is for if there is no group.of.interest
      # equivilant to group.of.interest = "all"
      # with nos?
      if (include.nos == FALSE) {
        which(subgroup == "NOS") -> nos.index
        subgroup <- subgroup[-nos.index]
        gene.exp <- gene.exp[-nos.index]
        nmb.name <- nmb.name[-nos.index]
        subgroup <- factor(subgroup,levels = subgroup.include)
        
        
      } else {
        subgroup <- factor(subgroup,levels = subgroup.include)
        
        
      }
      
      #sort out these titles
      
      if (is.character(hgnc.id[1,]) == "TRUE") {
        title = paste0(gene.of.interest," (",hgnc.id,") expression over subgroups")
      }
      else if (is.character(hgnc.id[1,]) == "FALSE") {
        title = paste0(gene.of.interest," expression over subgroups")
      }
      
    }
    
    
    # standard scale?
    if (standard.scale == TRUE) {
      scale <- c(-1,21)
    }else{
      scale <- c(min(gene.exp - 0.5),max(gene.exp) + 0.5)
    }
    
    
    ### numbers and colours sorted, hardcoded for now
    # how many of each subgroup?
    table(subgroup) -> numbers
    names = paste0(levels(subgroup), " (n=", numbers, ")")
    
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
    
    ###plotting happens
    
    par(mar = c(5, 5, 5, 6), las = 1)
    par(xpd = FALSE)
    # boxplot of expression vs subgroup
    boxplot(
      gene.exp ~ subgroup, col = sbgrp.col, ylim = scale, ylab = "Expression (VSD)", names = names, main = title
    )
    # ablind on the median
    abline(h = median(gene.exp), lty = 2)
    # plot outside the margin
    par(xpd = TRUE)
    # write median using the largest x value with an adjustment
    text(par("usr")[2],median(gene.exp),"Median", cex = 1, adj = -.1)
    
    
    ### ANOVA and text
    # so we need to check the number of subgroups to see what is appropriate, ANOVA or T TEST
    # more than 2? ANOVA
    if (length(unique(subgroup)) > 2) {
      # this cat is for the log file
      cat("\nANOVA of ",  gsub(paste0(gene.of.interest, " \\(|) expression"), "", title),":\n\n")
      
      # do the ANOVA and summarise it
      aov(gene.exp ~ subgroup) -> stat
      summary(stat) -> stat.summary
      
      # print for the log
      print(stat.summary)
      
      # pull out the values F/p
      stat.summary[[1]]$"F value"[1] -> F
      stat.summary[[1]]$"Pr(>F)"[1] -> p.val
      
      # tweak the p value if it is tin or round it
      if (p.val < 0.001 & !is.na(p.val)) {
        p.text <- "p < 0.001"
      } else {
        p.text <- paste("p =" , round(p.val,3))
      }
      
      #round the f value
      f.text <- paste("F =",round(F,3))
      
      # write the text as above but under the other values
      text(par("usr")[2],median(gene.exp) - ((par("usr")[4] - par("usr")[3]) *
                                               .05), f.text, adj = -.1)
      text(par("usr")[2],median(gene.exp) - ((par("usr")[4] - par("usr")[3]) *
                                               .1), p.text, adj = -.1)
      
      ### T test
      # if there are just two groups...T TEST
    } else if (length(unique(subgroup)) == 2) {
      cat("\nT Test of ",  gsub(paste0(gene.of.interest, " \\(|) expression"), "", title),":\n\n")
      
      
      
      # do the t test and print the summary for the log
      t.test(gene.exp ~ subgroup) -> stat
      print(summary(stat))
      
      
      # same as above pull out the df/p val and round appropriately
      stat[[2]] -> df
      stat[[3]] -> p.val
      
      
      if (p.val < 0.001 & !is.na(p.val)) {
        p.text <- "p < 0.001"
      } else {
        p.text <- paste("p =" , round(p.val,3))
      }
      f.text <- paste("DF =",signif(df,3))
      
      
      # write these on the graph
      text(par("usr")[2],median(gene.exp) - ((par("usr")[4] - par("usr")[3]) *
                                               .05), f.text, adj = -.1)
      text(par("usr")[2],median(gene.exp) - ((par("usr")[4] - par("usr")[3]) *
                                               .1), p.text, adj = -.1)
      
      
    }
    # this bit is for log file creating a data frame of what was plotted
    gene.exp.output <- rbind(gene.exp, as.character(subgroup))
    rownames(gene.exp.output) <-
      c(paste0(hgnc.id," (", gene.of.interest, ")"), "Subgroup")
    colnames(gene.exp.output) <- nmb.name
    #gene.exp.output
    # turn the sink off
    LogSinker(start = FALSE)
    # the output is returned for saving to a csv in the Rmd file
    return(gene.exp.output)
    
    
  }
