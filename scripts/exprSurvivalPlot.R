exprSurvivalPlot <- function(gene.of.interest, vsd, annot, subgroup.include = "all", include.nos = FALSE, hgnc.id = hgnc.id, flag.limit = 0.75, output.dirname , bin,  cont = FALSE) {
    
  #
  #This function creates KM survival plots for all subgroups or individual ones with accompanying statistics
  #
  # gene.of.interest = gene of interest ENSEMBL
  # vsd = expression vsd file
  # annot= expression annotation file
  # subgroup.include =  which subgroups to include but used to tell which data is INCLUDED in the plot not for comparison, all, WNT, SHH, Grp3, Grp4
  #^^^^^ case sensitive! ^^^^^
  # include.nos = include NOS samples (TRUE)
  # hgnc.id = the gene id created earlier
  # flag limit = the limit to give a warning as more than that percentage are not expressed
  #
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
  which(gene.of.interest == gene.names) -> gene.match

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
    if(bin == "median"){
    # get the median vsd for the gene of interest
    median(gene.exp[keep.index]) -> median.gene.of.interest.vsd
    bin.name <- paste0(" Binned by Median (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")     
    } else if(bin == "mean"){
    mean(gene.exp[keep.index]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by Mean (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    } else if(bin == "quart.100"){
     unname(quantile(gene.exp[keep.index])[5]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 100% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    } else if(bin == "quart.75") {
      unname(quantile(gene.exp[keep.index])[4]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 75% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    } else if(bin == "quart.50") {
      unname(quantile(gene.exp[keep.index])[3]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 50% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    } else if(bin == "quart.25") {
      unname(quantile(gene.exp[keep.index])[2]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 25% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")   
    } else if(bin == "quart.0") {
      unname(quantile(gene.exp[keep.index])[1]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 0% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")   
    } else if(grepl("split", bin) == TRUE){
      bin.name <- paste0(" Binned by <25%, 25-75% and >75% (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")  
    } else {
      max(gene.exp[keep.index])*bin-> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by ", bin, "% (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    }
      

      
     if(grepl("split", bin)) {
       
      low <- strsplit(bin, split = "_")[[1]][2]
      high <- strsplit(bin, split = "_")[[1]][3]
      
      quantiles <- quantile(gene.exp[keep.index])
    
      quantiles[grep(low, names(quantiles))]
      quantiles[grep(high, names(quantiles))]
      
      unname(quantile(gene.exp[keep.index])[2])
      
      gene.exp.bin <- ifelse(gene.exp[keep.index] <  quantiles[grep(low, names(quantiles))],
             "Low", 
             ifelse(gene.exp[keep.index] > quantiles[grep(low, names(quantiles))] & quantiles[grep(high, names(quantiles))] < gene.exp[keep.index], 
                    "Middle", "High"))
      gene.exp.bin <- factor(gene.exp.bin, levels = c("Low","Middle", "High"))

      bin.name <- paste0(" Binned by <25%, 25-75% and >75% (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")     
       
       
     } else {
    # use this to assign either high or low, then make it into a factor
    ifelse(gene.exp[keep.index] > median.gene.of.interest.vsd,"High","Low") -> gene.exp.bin
    factor(gene.exp.bin, levels = c("Low","High")) -> gene.exp.bin
     }
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
      title <- paste("Survival for", hgnc.id[1,1], bin.name)
    } else {
      title <-
        paste0("Survival for ", hgnc.id[1,1], " (", paste0(temp.subgroup.include, collapse = ", "),")",bin.name)
    }
    
    # create the Kaplan-Meier
    KM <-      try(survfit(Surv(time.os[keep.index],status.os[keep.index]) ~ gene.exp.bin, type = "kaplan-meier", conf.type = "log"))
    
    # output flag and KM data to log file
    cat("Flag Limit: ", flag.limit, "\n\n")
    cat("\nKM data for: ", hgnc.id[[1]], "\n\n")
    print(KM)
    
    # plot the kaplan-meier plot
    # hardcoded red blue
    try(plot(
      KM, xlab = "Time (years)", ylab = "Cumulative Overall Survival", lwd = 2, col =
        km.colours, main = title, mark.time= TRUE 
    ))
    # add the legend
    
    
    legend(
      x = "topright", col = km.colours, lwd = 2, legend = levels(as.factor(gene.exp.bin)), bg = "white"
    )
    # do a Log-Rank test
    try(survdiff(Surv(time.os[keep.index],status.os[keep.index])  ~ gene.exp.bin)) -> surv.log.rank
    
    # output log rank
    cat("\nLog Rank data for: ", hgnc.id[[1]], "\n\n")
      print(surv.log.rank)
        
    # add the rounded p-value to the plot
      
    if(!is(surv.log.rank, "try-error")){
    1 - pchisq(surv.log.rank$chisq, length(surv.log.rank$obs) - 1) -> surv.p.val
    
    # output the p value
    cat("\nRounded P Value for: ", hgnc.id[[1]], "\n\n")
    cat(surv.p.val)
    
    
    if (surv.p.val < 0.001 & !is.na(surv.p.val)) {
      p.text <- "p < 0.001"
    } else {
      p.text <- paste(" p =" , round(surv.p.val,3))
    }
    
    }
    
    # create the numbers for the plot
    text(par("usr")[2],0.1,paste("N = ", sum(KM$n), p.text), adj = 1.1)
    # turn the sink off and return the output
    LogSinker(start = FALSE)
    return(surv.output)
  }


exprCoxPlot <-  function(gene.of.interest, vsd, annot, subgroup.include = "all", include.nos = FALSE, hgnc.id = hgnc.id, flag.limit = 0.75, bin, output.dirname = output.dirname, cont = FALSE) {
   
  #
  #This function creates Cox models for all subgroups or individual ones with accompanying statistics
  #
  # gene.of.interest = gene of interest ENSEMBL
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
    which(gene.of.interest == gene.names) -> gene.match
    
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
    
    if(bin == "median"){
      # get the median vsd for the gene of interest
      median(gene.exp[keep.index]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by Median (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")     
    } else if(bin == "mean"){
      mean(gene.exp[keep.index]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by Mean (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    } else if(bin == "quart.100"){
      unname(quantile(gene.exp[keep.index])[5]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 100% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    } else if(bin == "quart.75") {
      unname(quantile(gene.exp[keep.index])[4]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 75% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    } else if(bin == "quart.50") {
      unname(quantile(gene.exp[keep.index])[3]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 50% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    } else if(bin == "quart.25") {
      unname(quantile(gene.exp[keep.index])[2]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 25% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")   
    } else if(bin == "quart.0") {
      unname(quantile(gene.exp[keep.index])[1]) -> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by 0% Quartile (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")   
    } else if(grepl("split", bin) == TRUE){
      bin.name <- paste0(" Binned by <25%, 25-75% and >75% (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")  
    } else {
      max(gene.exp[keep.index])*bin-> median.gene.of.interest.vsd
      bin.name <- paste0(" Binned by ", bin, "% (", sprintf("%.2f", round(median.gene.of.interest.vsd,2)), ")")      
    }
    
    if(grepl("split", bin)) {
      
      low <- strsplit(bin, split = "_")[[1]][2]
      high <- strsplit(bin, split = "_")[[1]][3]
      
      quantiles <- quantile(gene.exp[keep.index])
      
      quantiles[grep(low, names(quantiles))]
      quantiles[grep(high, names(quantiles))]
      
      unname(quantile(gene.exp[keep.index])[2])
      
      gene.exp.bin <- ifelse(gene.exp[keep.index] <  quantiles[grep(low, names(quantiles))],
                             "Low", 
                             ifelse(gene.exp[keep.index] > quantiles[grep(low, names(quantiles))] & quantiles[grep(high, names(quantiles))] < gene.exp[keep.index], 
                                    "Middle", "High"))
      gene.exp.bin <- factor(gene.exp.bin, levels = c("Low","Middle", "High"))
      
         
      
      
    } else {
      # use this to assign either high or low, then make it into a factor
      ifelse(gene.exp[keep.index] > median.gene.of.interest.vsd,"High","Low") -> gene.exp.bin
      factor(gene.exp.bin, levels = c("Low","High")) -> gene.exp.bin
    }
    #length(which(gene.exp[keep.index] == 0)) / length(gene.exp[keep.index]) -> perc.non.express
    
    
    #### cox model categorical
    cat.cox.model <- try(coxph(Surv(time.os[keep.index],status.os[keep.index])  ~ gene.exp.bin))
    
    # print the output to log
    cat("\n\nCox model categorical: ", hgnc.id[[1]], "\n\n")
    print(cat.cox.model)
    cat.coef <- try(coef(cat.cox.model))
    # this output table will be used in the kable to create the plot and massage into the right format
    output.table <- data.frame(
    "n" = paste0(cat.cox.model$nevent, "/", cat.cox.model$n), 
    "HR" = round(exp(cat.coef), 3),
    "CI" = paste0(signif(summary(cat.cox.model)$conf.int[[3]] ,3), "-" ,signif(summary(cat.cox.model)$conf.int[[4]],3)),
    "P value" =  signif(summary(cat.cox.model)$coefficients[[5]], 3)
    )
    colnames(output.table)[3:4]<- as.character(c("95% CI", "P value"))
    

    rownames(output.table) <- paste0(gsub("gene.exp.bin", "", rownames(output.table)), " Expression (Cat.)")
    
    
    # library(xtable)
    #x.table <- xtable(output.table)
    #sink(paste0(output.dir,paste0("/",hgnc.id,".cox.model.cat.html")))
    #print(x.table,type="html")
    #sink()
    
    #### cox model continuous
    cont.cox.model <- try(coxph(Surv(time.os[keep.index],status.os[keep.index])  ~ gene.exp[keep.index]))
    cont.coef <- try(coef(cont.cox.model))
    
    # print for the log
    cat("\n\nCox model continuous: ", hgnc.id[[1]], "\n\n")
    print(cont.cox.model)
    
    #second half of the output table, massage as above
    output.table2 <- data.frame(
    "n" = paste0(cont.cox.model$nevent, "/", cont.cox.model$n), 
    "HR" = round(exp(cont.coef), 3),
    "CI" = paste0(signif(summary(cont.cox.model)$conf.int[[3]] ,3), "-" ,signif(summary(cont.cox.model)$conf.int[[4]],3)),
    "p value" =  signif(summary(cont.cox.model)$coefficients[[5]], 3)
    )
    colnames(output.table2)[3:4]<- as.character(c("95% CI", "P value"))
    
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
