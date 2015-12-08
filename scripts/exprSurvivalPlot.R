exprSurvivalPlot <- function(goi, vsd, annot, subgroup.include = "all", include.nos = FALSE, hgnc.id = hgnc.id, flag.limit = 0.75, output.dirname ,  cont = FALSE) {
    
  #
  #This function creates KM survival plots for all subgroups or individual ones with accompanying statistics
  #
  # goi = gene of interest ENSEMBL
  # vsd = expression vsd file
  # annot= expression annotation file
  # subgroup.include =  which subgroups to include but used to tell which data is INCLUDED in the plot not for comparison, all, WNT, SHH, Grp3, Grp4
  #^^^^^ case sensitive! ^^^^^
  # include.nos = include NOS samples (TRUE)
  # hgnc.id = the gene id created earlier
  # flag limit = the limit to give a warning as more than that percentage are not expressed
  #
  # output.dirname = location for the log file
  # cont = whether to continue the log file or make a new one (when plotting subgroups)
  #
  
  # create/append to the log file
  LogSinker(
      data.type = "Survival", outdir = paste0(output.dirname,"/Survival/"), start = TRUE, cont = cont
    )
    
  
  
  # trim off any .1 from the end of ensembl IDs
  gsub('\\.[0-9]+$', '', rownames(vsd)) -> gene.names
  # find the gene match from the data set and get the expression for that
  which(goi == gene.names) -> gene.match

   # stop if no match found 
    if (length(gene.match) < 1) {
      stop("no gene match")
    }else if (length(gene.match) > 2) {
      stop("more than one gene matches")
    }
    
  
  # to not include NOS 
    if (include.nos == FALSE) {
      
      # pull out the data and remove nos samples
      annot$Group -> subgroup
      which(subgroup == "NOS") -> nos.index
      subgroup <- subgroup[-nos.index]
      
      # pull out time to PFS, OS, STATUS etc.
      annot$time_pfs[-nos.index] -> time.pfs
      annot$time_folow_up[-nos.index] -> time.os
      annot$Status_follow_up[-nos.index] -> status.os
      annot$Relapse[-nos.index] -> status.pfs
      # create died of disease true false vector
      status.os <- status.os == "DOD"
      # Progression/Relapse TRUE/FALSE vector
      status.pfs <- status.pfs == "Progression" | status.pfs == "Relapse"
      # gene expression
      as.numeric(vsd[gene.match,-nos.index]) -> gene.exp
      Index <- annot$Index[-nos.index]
      # pull out NMB names
      
      # sort out the subgroups to include
      if (paste0(subgroup.include, collapse = "|") == "all") {
        subgroup.data -> temp.subgroup.include
      } else {
        subgroup.include -> temp.subgroup.include
        
      }
    } else {
      
      # same as above but without removing the NOS values
      annot$Group -> subgroup
      annot$time_pfs -> time.pfs
      annot$time_folow_up -> time.os
      annot$Status_follow_up -> status.os
      status.os <- status.os == "DOD"
      annot$Relapse -> status.pfs
      status.pfs <- status.pfs == "Progression" | status.pfs == "Relapse"
      as.numeric(vsd[gene.match,]) -> gene.exp
      Index <- annot$Index
      
      if (paste0(subgroup.include, collapse = "|") == "all") {
        c(subgroup.data,"NOS") -> temp.subgroup.include
      } else {
        subgroup.include -> temp.subgroup.include
        
      }
      
    }
    
    
    
    # for the log file 
    cat("\nSubgroups: ", temp.subgroup.include, "\n\n")
    
    # creaste the subgroup keep index
    which(subgroup %in% temp.subgroup.include) -> keep.index
    
    
    ### divide into two ? is this alway appropriate tertiles?
    # currently hard coded to split by median but this will be changing 
    
    # get the median vsd for the gene of interest
    median(gene.exp[keep.index]) -> median.goi.vsd
    # use this to assign either high or low, then make it into a factor
    ifelse(gene.exp[keep.index] > median.goi.vsd,"high","low") -> gene.exp.bin
    factor(gene.exp.bin, levels = c("low","high")) -> gene.exp.bin
    
    # this is the output for records, creates a data frame
    surv.output <- data.frame(
      PFS = time.pfs[keep.index],
      OS = time.os[keep.index],
      Status = status.os[keep.index],
      Subgroup = subgroup[keep.index],
      Bin = gene.exp.bin,
      Exprs = gene.exp[keep.index]
    )
    
    # transpose the output and rename. just for aesthetics
    surv.output <- t(surv.output)
    colnames(surv.output) <- Index[keep.index]
    
    
    # find out how many samples do not express, if greater than the flag value you get  a warning message
    length(which(gene.exp[keep.index] == 0)) / length(gene.exp[keep.index]) -> perc.non.express
    
    if (perc.non.express > flag.limit) {
      non.exprs.flag = TRUE
    }else{
      non.exprs.flag = FALSE
    }
    
    if (non.exprs.flag) {
      text(3,0.3, paste(
        round(perc.non.express * 100,1), "% of sample \n do not express!"
      ), pos = 4)
    }
    
    # sort out the names
    if (paste0(subgroup.include, collapse = "|") == "all") {
      title <- paste("Survival for", hgnc.id)
    } else {
      title <-
        paste0("Survival for ", hgnc.id, " (", paste0(temp.subgroup.include, collapse = ", "),")")
    }
    
    # create the Kaplan-Meier
    KM <-
      survfit(Surv(time.os[keep.index],status.os[keep.index]) ~ gene.exp.bin, type = "kaplan-meier", conf.type = "log")
    
    # output flag and KM data to log file
    cat("Flag Limit: ", flag.limit, "\n\n")
    cat("\nKM data for: ", hgnc.id[[1]], "\n\n")
    print(KM)
    
    # plot the kaplan-meier plot
    # hardcoded red blue
    plot(
      KM, xlab = "Time (years)", ylab = "Cumulative Overall Survival", lwd = 2, col =
        km.colours, main = title
    )
    # add the legend
    
    legend(
      x = "topright", col = km.colours, lwd = 2, legend = levels(as.factor(gene.exp.bin)), bg = "white"
    )
    # do a Log-Rank test
    survdiff(Surv(time.os[keep.index],status.os[keep.index])  ~ gene.exp.bin) -> surv.log.rank
    
    # output log rank
    cat("\nLog Rank data for: ", hgnc.id[[1]], "\n\n")
        print(surv.log.rank)
        
    # add the rounded p-value to the plot
    1 - pchisq(surv.log.rank$chisq, length(surv.log.rank$obs) - 1) -> surv.p.val
    
    # output the p value
    cat("\nRounded P Value for: ", hgnc.id[[1]], "\n\n")
    cat(surv.p.val)
    
    
    # create the numbers for the plot
    text(par("usr")[2],0.1,paste("N = ", sum(KM$n)," p =",round(surv.p.val, 3)), adj = 1.1)
    # turn the sink off and return the output
    LogSinker(start = FALSE)
    return(surv.output)
  }


exprCoxPlot <-  function(goi, vsd, annot, subgroup.include = "all", include.nos = FALSE, hgnc.id = hgnc.id, flag.limit = 0.75, output.dirname = output.dirname, cont = FALSE) {
   
  #
  #This function creates Cox models for all subgroups or individual ones with accompanying statistics
  #
  # goi = gene of interest ENSEMBL
  # vsd = expression vsd file
  # annot= expression annotation file
  # subgroup.include =  which subgroups to include but used to tell which data is INCLUDED in the plot not for comparison, all, WNT, SHH, Grp3, Grp4
  #^^^^^ case sensitive! ^^^^^
  # include.nos = include NOS samples (TRUE)
  # hgnc.id = the gene id created earlier
  # flag limit = the limit to give a warning as more than that percentage are not expressed
  #
  # output.dirname = location for the log file
  # cont = whether to continue the log file or make a new one (when plotting subgroups)
  #
  # This function is sperate to make the knitting easier in the report it could be done easily above 
  #
  # turn on the sink, with option to append, this would be appended
  
   LogSinker(
      data.type = "Survival", outdir = paste0(output.dirname,"/Survival/"), start = TRUE, cont = cont
    )
    
    # trim off any .1 from the end of ensembl IDs
    gsub('\\.[0-9]+$', '', rownames(vsd)) -> gene.names
    # find the gene match from the data set and get the expression for that
    which(goi == gene.names) -> gene.match
    
    # stop if no match found 
    if (length(gene.match) < 1) {
      stop("no gene match")
    }else if (length(gene.match) > 2) {
      stop("more than one gene matches")
    }
    
    
    # to not include NOS 
    if (include.nos == FALSE) {
      
      # pull out the data and remove nos samples
      annot$Group -> subgroup
      which(subgroup == "NOS") -> nos.index
      subgroup <- subgroup[-nos.index]
      
      # pull out time to PFS, OS, STATUS etc.
      annot$time_pfs[-nos.index] -> time.pfs
      annot$time_folow_up[-nos.index] -> time.os
      annot$Status_follow_up[-nos.index] -> status.os
      annot$Relapse[-nos.index] -> status.pfs
      # create died of disease true false vector
      status.os <- status.os == "DOD"
      # Progression/Relapse TRUE/FALSE vector
      status.pfs <- status.pfs == "Progression" | status.pfs == "Relapse"
      # gene expression
      as.numeric(vsd[gene.match,-nos.index]) -> gene.exp
      Index <- annot$Index[-nos.index]
      # pull out NMB names
      
      # sort out the subgroups to include
      if (paste0(subgroup.include, collapse = "|") == "all") {
        subgroup.data -> temp.subgroup.include
      } else {
        subgroup.include -> temp.subgroup.include
        
      }
    } else {
      
      # same as above but without removing the NOS values
      annot$Group -> subgroup
      annot$time_pfs -> time.pfs
      annot$time_folow_up -> time.os
      annot$Status_follow_up -> status.os
      status.os <- status.os == "DOD"
      annot$Relapse -> status.pfs
      status.pfs <- status.pfs == "Progression" | status.pfs == "Relapse"
      as.numeric(vsd[gene.match,]) -> gene.exp
      Index <- annot$Index
      
      if (paste0(subgroup.include, collapse = "|") == "all") {
        c(subgroup.data,"NOS") -> temp.subgroup.include
      } else {
        subgroup.include -> temp.subgroup.include
        
      }
      
    }
    
    which(subgroup %in% temp.subgroup.include) -> keep.index
    
    
    # creaste the subgroup keep index
    which(subgroup %in% temp.subgroup.include) -> keep.index
    
    
    ### divide into two ? is this alway appropriate tertiles?
    # currently hard coded to split by median but this will be changing 
    
    # get the median vsd for the gene of interest
    median(gene.exp[keep.index]) -> median.goi.vsd
    # use this to assign either high or low, then make it into a factor
    ifelse(gene.exp[keep.index] > median.goi.vsd,"high","low") -> gene.exp.bin
    factor(gene.exp.bin, levels = c("low","high")) -> gene.exp.bin
    
    #length(which(gene.exp[keep.index] == 0)) / length(gene.exp[keep.index]) -> perc.non.express
    
    
    #### cox model categorical
    cox.model <- coxph(Surv(time.os[keep.index],status.os[keep.index])  ~ gene.exp.bin)
    
    # print the output to log
    cat("\n\nCox model categorical: ", hgnc.id[[1]], "\n\n")
    print(cox.model)
    coef <- coef(cox.model)
    # this output table will be used in the kable to create the plot and massage into the right format
    output.table <-
      data.frame(
        coef = coef, exp.coef = exp(coef), se.coef = summary(cox.model)$coefficients[[3]], z = summary(cox.model)$coefficients[[4]], p.val = summary(cox.model)$coefficients[[5]]
      )
    row.names(output.table) <- "High Expression (Cat.)"
    
    
    # library(xtable)
    #x.table <- xtable(output.table)
    #sink(paste0(output.dir,paste0("/",hgnc.id,".cox.model.cat.html")))
    #print(x.table,type="html")
    #sink()
    
    #### cox model continuous
    cox.model <- coxph(Surv(time.os[keep.index],status.os[keep.index])  ~ gene.exp[keep.index])
    
    # print for the log
    cat("\n\nCox model continuous: ", hgnc.id[[1]], "\n\n")
    print(cox.model)
    
    #second half of the output table, massage as above
    output.table2 <-
      data.frame(
        coef = coef, exp.coef = exp(coef), se.coef = summary(cox.model)$coefficients[[3]], z = summary(cox.model)$coefficients[[4]], p.val = summary(cox.model)$coefficients[[5]]
      )
    row.names(output.table2) <- "High Expression (Cont.)"
    
    large.table <- rbind(output.table,output.table2)
    
    # return the table
    print(large.table)
    #x.table <- xtable(output.table)
    #sink(paste0(output.dir,paste0("/",hgnc.id,".cox.model.cont.html")))
    # print(x.table,type="html")
    #sink()
    
    
    # turn the sink off
    LogSinker(start = FALSE)
    
    # return the table
    return(large.table)
  }
