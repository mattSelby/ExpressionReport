

#
# leaving this bit of pre processing at the start of the script
# I may make it into a function later but I'm not sure yet
#
# it sorts the data out and then pulls out the gene.of.interest specific data for all the variant bits
# reordered but not sure if that is necessary
#

###        do we need this or can we just use the annot$index?

### removes the .XX trailing ensembl.id and extracts simple sample name
gsub("../output/Sample_","",cuff.samp$file) -> sample.names
gsub("_sortedRG.cxb","",sample.names) -> sample.names
gsub("_","",sample.names) -> sample.names
match(annot$Index,sample.names) -> cufflinks.reordering.index
gsub('\\.[0-9]+$', '', isoforms.attr$nearest_ref_id) -> gene.names

colnames(isoforms.counts) <- sample.names
colnames(genes.counts) <- sample.names
colnames(cds.counts) <- sample.names
colnames(tss.counts) <- sample.names

colnames(isoforms.fpkm) <- sample.names
colnames(genes.fpkm) <- sample.names
colnames(cds.fpkm) <- sample.names
colnames(tss.fpkm) <- sample.names



### removed the correction step
#match(annot$Index,sample.names) -> cufflinks.reordering.index

#round(t(apply(isoforms.counts, 1, function(x){x*(cuff.samp$internal_scale)})),0) -> isoforms.count.corrected
#round(t(apply(genes.counts, 1, function(x){x*(cuff.samp$internal_scale)})),0) -> genes.count.corrected
#round(t(apply(tss.counts, 1, function(x){x*(cuff.samp$internal_scale)})),0) -> tss.count.corrected
#round(t(apply(cds.counts, 1, function(x){x*(cuff.samp$internal_scale)})),0) -> cds.count.corrected

gsub('\\.[0-9]+$', '', isoforms.attr$nearest_ref_id) -> ensembl.ids
### put a unique in to stop the error message
unique(as.character(isoforms.attr[which(ensembl.ids == gene.of.interest),"gene_id"])) -> gene.of.interest.id


#### find the gene/isoform/cds/tss
genes.attr[which(genes.attr$gene_id == gene.of.interest.id),"tracking_id"] -> genes.index
cds.attr[which(cds.attr$gene_id == gene.of.interest.id),"tracking_id"] -> cds.index
isoforms.attr[which(isoforms.attr$gene_id == gene.of.interest.id),"tracking_id"] -> isoforms.index
tss.attr[which(tss.attr$gene_id == gene.of.interest.id),"tracking_id"] -> tss.index

### pull out the data
gene.of.interest.gene.fpkm <- genes.fpkm[genes.index,cufflinks.reordering.index]
gene.of.interest.isoform.fpkm <- isoforms.fpkm[isoforms.index,cufflinks.reordering.index]
gene.of.interest.isoform.attr <- isoforms.attr[isoforms.index,]
gene.of.interest.cds.fpkm <- cds.fpkm[cds.index,cufflinks.reordering.index]
gene.of.interest.tss.fpkm <- tss.fpkm[tss.index,cufflinks.reordering.index]

gene.of.interest.gene.counts <- genes.counts[genes.index,cufflinks.reordering.index]
gene.of.interest.isoform.counts <- isoforms.counts[isoforms.index,cufflinks.reordering.index]
gene.of.interest.cds.counts <- cds.counts[cds.index,cufflinks.reordering.index]
gene.of.interest.tss.counts <- tss.counts[tss.index,cufflinks.reordering.index]


#### subgroups and also get rid of NOS


### function to assign colours


assignColors <- function(fact) {
  
  #
  # fact = the thing to make colours for, important bit is the length
  # function to create colours for the plotting here
  # if there are more than 8 colours (the extent the package goes to) then repeat the colours
  #
  
  # unique values of the input
  n <- length(unique(fact))
  if (n > 8) {
    rep(brewer_pal(palette = "Pastel1")(8),10)[1:n] -> brewer.cols
  }else{
    brewer_pal(palette = "Pastel1")(n) -> brewer.cols
  }
  # create a dataframe of the colours and add grey for NA
  fact.col <-
    data.frame(level.name = unique(fact),level.colors = brewer.cols)
  as.character(fact.col$level.colors) -> out.cols
  out.cols[is.na(fact.col$level.name)] <- "grey"
  out.cols[match(fact,fact.col$level.name)] -> out
  return(out)
}



sortVariantData <- function(subgroup.include = "all", gene.of.interest.isoform.fpkm, gene.of.interest.cds.fpkm ,gene.of.interest.tss.fpkm, 
                              gene.of.interest.isoform.attr, gene.of.interest.cds.attr, gene.of.interest.tss.attr, include.nos = FALSE, output.dirname = output.dirname) {
     
  #
  # This function is to sort out the transcript data, it will return a list of 
  # lists with the processed data has been split for ease of plotting etc. in
  # the .Rmd document
  #
  # subgroup.include =  which subgroups to include but used to tell which data
  # is INCLUDED in the plot not for comparison, all, WNT, SHH, Grp3, Grp4
  #^^^^^ case sensitive! ^^^^^
  # THIS WILL AFFECT ALL OF THE LATER PLOTS
  # gene.of.interest.isoform.fpkm = fpkm data for isoforms
  # gene.of.interest.cds.fpkm = cds data for isoforms
  # gene.of.interest.tss.fpkm = tss data for isoforms
  # gene.of.interest.isoform.attr = attributes of the isoform data
  # gene.of.interest.cds.attr = attributes of the cds data
  # gene.of.interest.tss.attr = attributes of the tss data
  # include.nos = include NOS samples (TRUE)
  # output.dirname = location for the log file
  #
  #
  
  # starts the log file
  
  LogSinker(
    data.type = "TranscriptProduction", outdir = paste0(output.dirname,"/Transcripts/"), start = TRUE, cont = FALSE
  )
  
  # assigns subgroup and sample name data
    annot$Group -> subgroup
    annot$Index -> sample.name
    
    # this bit will remove any nos samples and prune all the relveant data
    if (include.nos == FALSE) {
      which(subgroup == "NOS") -> nos.index
      subgroup <- subgroup[-nos.index]
      gene.of.interest.isoform.fpkm <- gene.of.interest.isoform.fpkm[-nos.index]
      gene.of.interest.cds.fpkm <- gene.of.interest.cds.fpkm[-nos.index]
      gene.of.interest.tss.fpkm <- gene.of.interest.tss.fpkm[-nos.index]
      gene.of.interest.isoform.attr <- gene.of.interest.isoform.attr[-nos.index]
      
    }
    
    # here we use the subgroup include to remove any subgroups which are not to be included in the analysis
    # collapse and search subgroup include for "all" then invert
    if (!paste0(subgroup.include, collapse = "|") == "all") {
      if (include.nos == TRUE) {
        #this adds an NOS if they are to be included
        subgroup.include <- c(subgroup.include, "NOS")
      }
      # which of the subgroups/NOS are to be kept?
      which(subgroup %in% subgroup.include) -> keep.index
      # remove the data then re level to correct plotting order
      subgroup <- subgroup[keep.index]
      subgroup <- factor(subgroup, levels = subgroup.include)
      
      
    } else if (include.nos == TRUE) {
      # this is for include NOS and all subgroups
      subgroup <-
        factor(subgroup,levels = c(subgroup.data, "NOS"));1:length(subgroup) -> keep.index
    } else {
      # this is for all subgroups
      subgroup <-
        factor(subgroup,levels = subgroup.data); 1:length(subgroup) -> keep.index
    }
    
    
    ### generate the isoform data
    # if to check that there is data present before proceeding
    if (!nrow(gene.of.interest.isoform.fpkm) == 0) {
    #create output list
      data.Isoforms <- list()
      # log2 +1 the data in a for loop putting it into the output list
      for (i in 1:nrow(gene.of.interest.isoform.fpkm)) {
        log2(as.numeric(gene.of.interest.isoform.fpkm[i,keep.index]) + 1) -> data.Isoforms[[i]]
      }
      # bind that together
      do.call(rbind, data.Isoforms) -> data.Isoforms
      
      ### find only isoform of interest
      
      #move the column names over
      colnames(gene.of.interest.isoform.fpkm[keep.index]) -> colnames(data.Isoforms)
      # find your gene.of.interest
      match(isoforms.attr$tracking_id, rownames(gene.of.interest.isoform.fpkm)) -> index
      # strip any na values out
      which(!is.na(index)) -> index
      # sort the ensembl name out
      gsub('\\.[0-9]+$', '',isoforms.attr[index,"nearest_ref_id"]) -> rownames(data.Isoforms)
      ### sort out names
      # take out the ENSG value as only interested in ENST
      grep("ENSG",rownames(data.Isoforms)) -> remove
      # sort the data by removing any that are not needed
      data.Isoforms <- data.Isoforms[-remove,,drop = FALSE]
    } else {
      # this shows none are found.....
      cat("\n\nISOFORM NOT FOUND\n")
      
    }
    
    
    ########### generate the cds data
    # if to check that there is data present before proceeding
    if (!nrow(gene.of.interest.cds.fpkm) == 0) {
      #create output list
      data.CDS <- list()
      # log2 +1 the data in a for loop putting it into the output list
      for (i in 1:nrow(gene.of.interest.cds.fpkm)) {
        log2(as.numeric(gene.of.interest.cds.fpkm[i,keep.index]) + 1) -> data.CDS[[i]]
      }
      # bind that together
      do.call(rbind, data.CDS) -> data.CDS
      
      
      ### sort names out
      #move the column names over
      colnames(gene.of.interest.cds.fpkm[keep.index]) -> colnames(data.CDS)
      # find your gene.of.interest
      match(cds.attr$tracking_id, rownames(gene.of.interest.cds.fpkm)) -> index
      # strip any na values out
      which(!is.na(index)) -> index
      # assign names to the data
      cds.attr[index,"tracking_id"] -> rownames(data.CDS)
    } else {
      # this shows none are found.....
      cat("\n\nCDS NOT FOUND\n")
      
    }
    
    ########### generate the tss data
    # if to check that there is data present before proceeding
    if (!nrow(gene.of.interest.tss.fpkm) == 0) {
      #create output list
      data.TSS <- list()
      # log2 +1 the data in a for loop putting it into the output list
      for (i in 1:nrow(gene.of.interest.tss.fpkm)) {
        log2(as.numeric(gene.of.interest.tss.fpkm[i,keep.index]) + 1) -> data.TSS[[i]]
      }
      # bind that together
      do.call(rbind, data.TSS) -> data.TSS
      
      ### sort names out
      #move the column names over
      colnames(gene.of.interest.tss.fpkm[keep.index]) -> colnames(data.TSS)
      # find your gene.of.interest
      match(tss.attr$tracking_id, rownames(gene.of.interest.tss.fpkm)) -> index
      # strip any na values out
      which(!is.na(index)) -> index
      # assign names to the data
      tss.attr[index,"tracking_id"] -> rownames(data.TSS)
      
    }  else {
      # this shows none are found.....
      cat("\n\nTSS NOT FOUND\n")
      
    }
    
    # this is to create your final outpt from this function
    output <- list()
    
    # check for isoform data, if yes put it into slot 1 of the list, if no put NA in
    if (!nrow(gene.of.interest.isoform.fpkm) == 0) {
      output[[1]] <- data.Isoforms
      
    } else {
      output[[1]] <- NA
      
    }
    
    # check for tss data, if yes put it into slot 2 of the list, if no put NA in
    if (!nrow(gene.of.interest.tss.fpkm) == 0) {
      output[[2]] <- data.TSS
      
    } else {
      output[[2]] <- NA
      
    }
    
    
    # check for cds data, if yes put it into slot 3 of the list, if no put NA in
    if (!nrow(gene.of.interest.cds.fpkm) == 0) {
      output[[3]] <- data.CDS
      
    } else {
      output[[3]] <- NA
      
    }
    
    # add the keep index to the list, slot 4
    output[[4]] <- keep.index
    # add your sample names to the list (- keep index), slot 5
    output[[5]] <- sample.name[keep.index]
    # add your subgroup to slot 6
    output[[6]] <- subgroup
    
    # turn off the sink
    LogSinker(
       start =  FALSE
    )
    
    return(output)
  }

### This is the sample data which is output and then used in further scripts....
### 
### slot [[1]] Isoform Data 
### slot [[2]] TSS Data 
### slot [[3]] CDS Data 
### slot [[4]] Keep Index Data
### slot [[5]] Sample Names Data 
### slot [[6]] Subgroup
###

###data<-output_for_isoforms[[1]]
###sample.names<-output_for_isoforms[[5]]
###subgroup<-output_for_isoforms[[6]]
###keep.index<-output_for_isoforms[[4]]

# 
# All of the plotting functions below use the output data from
# sort.variant.data() the different slots are input then used to plot each bit
# 

plotData <-function(data, sample.names, subgroup, keep.index, standard.scale = TRUE, multi = FALSE, group.of.interest = "all",output.dirname = output.dirname) {
  
  # 
  # not the most descriptive name but this is used to do any of the plotting
  # using variant data to create boxplots, no sink is used in this function as
  # the output is captured in sort.variant.data()
  # 
  # data = data to be plotted, e.g. isoform, cds
  # sample.names = slot [[5]] Sample Names Data
  # subgroup = slot [[6]] Subgroup
  # keep.index = slot [[4]] Keep Index Data
  # standard.scale = TRUE/FALSE (do we want a standard scale, this falls downa 
  # bit here as this is FPKM rather than vsd, may or may not need it??)
  # multi = do we want to plot at the bars on one chart (will be used as the
  # first plot) or do we want to split them all individually
  # group.of.interest = "all", which subgroups to compare for the stats, all or e.g. WNT vs others
  # output.dirname = output directory used for the sink in multi portion of the plotter
  #
  # pretty sure there stat bit of this could be functioned.....need to do that
  # the way this plots is set up for use with the Rmd
  
  # this identifies what sort of data we are dealing with and creates a
  # data.type for that
  
      if (any(grep("^ENST", rownames(data)))) {
      data.type = "Isoform"
    } else if (any(grep("^P", rownames(data)))) {
      data.type = "CDS"
    }else if (any(grep("^TSS", rownames(data)))) {
      data.type = "TSS"
    }
    
    
    # sort out the subgroup data frame for plotting
    temp.subgroup <- as.data.frame(subgroup)
    rownames(temp.subgroup) <- sample.names
    
    #### stick together the gene data which is plotted first and then the data of interest
    # ENSEMBL gene always goes at the start of plots regardless of osoform/cds/tss so it is added here
    # log 2 corrected
    # ensembl is put into the first row
    # data is then melted to make it plotable
    rbind(log2(as.numeric(gene.of.interest.gene.fpkm) + 1)[keep.index],data) -> temp.data
    rownames(temp.data)[1] <- gene.of.interest
    reshape2::melt(temp.data) -> temp
    
    # sort out an .1 which are ofund
    levels(temp$Var1) <- gsub('\\.[0-9]+$', '',levels(temp$Var1))
    
    
    ### sort out which subgroups samples are
    subgroup[match(temp$Var2,rownames(temp.subgroup))] -> multi.subgroup
    # re order subgroup so they will plot properly
    temp <- cbind(temp,Subgroup = factor(as.character(multi.subgroup),levels = levels(subgroup)))
    
    ### Group of interest sorted here, with title
    
    # is the group.of.interest anything other than all?
    if (!any(group.of.interest == "all")) {
      # add the groups and search to see which subgroup group.of.interest is
      groups <- c("WNT", "SHH", "Grp3", "Grp4")
      groups <-
        groups[grep(paste0(group.of.interest, collapse = "|"), groups, invert = TRUE)]
      # add pipes inbetween for later searching
      groups <- paste0(groups, collapse = "|")
      # are there NOS in the data, if so add that in
      if (any(levels(temp$Subgroup) == "NOS")) {
        final.level <-
          c("Other", group.of.interest, "NOS")
      } else {
        # or just add other to the other levels
        final.level <- c("Other", group.of.interest)
      }
      # search any groups that are need to be changed to other, e.g. if group.of.interest= 
      # Grp3 all others turne to Other
      temp$Subgroup <-
        factor(gsub(groups,"Other",temp$Subgroup), levels = final.level)
      
      # is there an hgnc.id? then make the title
      # is there a group.of.interest?
      if (is.character(hgnc.id[1,]) == "TRUE") {
        title = paste0(hgnc.id," ", data.type, " expression (",paste(group.of.interest, collapse  =
                                                                       ", "),") vs Others")
      }
      else if (is.character(hgnc.id[1,]) == "FALSE") {
        title = paste0(gene.of.interest," ",data.type,  " expression (",paste(group.of.interest, collapse  =
                                                                   ", "),") vs Others")
      }
      
      # if there is no group.of.interest do this
    } else {
      if (is.character(hgnc.id[1,]) == "TRUE") {
        title = paste0(hgnc.id," ", data.type, " expression over subgroups")
      }
      else if (is.character(hgnc.id[1,]) == "FALSE") {
        title = paste0(gene.of.interest," ",data.type," expression over subgroups")
      }
      
    }
    
    ### position parts for the plotter: for small/big labels and lines
    #
    # these are just positions to be used in plotting the various axis/ticks
    #
    # positions are the positions for writing the subgroup (minus the gap positions)
    positions <-
      1:(length(levels(temp$Var1)) * (length(levels(temp$Subgroup)) + 1))
    remove.positions <-
      seq(length(levels(temp$Subgroup)) + 1,max(positions),length(levels(temp$Subgroup)) +
            1)
    positions[-remove.positions] -> positions
    
    # this is for plotting the variant name below the subgroup
    
    mid.positions <-
      seq((length(levels(temp$Subgroup)) + 1) / 2,length(levels(temp$Var1)) *
            (length(levels(temp$Subgroup)) + 1),(length(levels(temp$Subgroup)) + 1))
    
    #seg.positions<-list()
    #for(i in 1:length(levels(temp$Var1))){
    #
    #  seg.positions[[i]]<-c(.5+((i-1)*(length(levels(temp$Subgroup))+1)),(length(levels(temp$Subgroup))+.5)+((i-1)*(length(levels(temp$Subgroup))+1)))
    #
    #
    #  }
    
    ### colours are sorted
    #
    # hardcoded but can be changed
    
    sbgrp.col <- c()
    if (any(levels(temp$Subgroup) == "Other")) {
      sbgrp.col  <- "grey"
    }
    if (any(levels(temp$Subgroup) == "WNT")) {
      sbgrp.col  <- c(sbgrp.col,"steelblue2")
    }
    if (any(levels(temp$Subgroup) == "SHH")) {
      sbgrp.col <- c(sbgrp.col,"tomato3")
    }
    if (any(levels(temp$Subgroup) == "Grp3")) {
      sbgrp.col <- c(sbgrp.col,"gold1")
    }
    if (any(levels(temp$Subgroup) == "Grp4")) {
      sbgrp.col <- c(sbgrp.col,"darkolivegreen1")
    }
    if (any(levels(temp$Subgroup) == "NOS")) {
      sbgrp.col <- c(sbgrp.col,"grey")
    }
    
    # standard scale? see above?
    if (standard.scale == TRUE) {
      scale <- c(-1,21)
    }else{
      scale <- c(min(temp$value - 0.5),max(temp$value) + 0.5)
    }
    
    ####
    ####
    ####  INDIVIDUAL PLOTTER HERE
    ####
    ####
    
    # this is for plotting all the data on one chart
    if (multi == FALSE) {
      # sort margin, xpd is here for debugging when plotting multiple times when checking to turn margin plotting off
      par(mar = c(13,5,8,2), xpd = FALSE)
      # boxplot of the log2 sample vs subgroup and gene/variant id
      boxplot(
        temp$value ~  temp$Subgroup + temp$Var1 ,
        col = sbgrp.col, ylim = scale,
        at = positions, cex.axis = 1.5, cex.main = 2, xaxt = 'n', las = 2, ylab = "log2(FPKM+1)", main = title
      )
      # ablines at intervals
      abline(h = seq(0,round(max(scale) / 2) * 2,2), lty = 2)
      
      
      # turn margin plotting on
      par(xpd = TRUE)
      ### sort out text sizes
      # this is done by the relative levels of the plotting object
      #
      # the more to plot the smaller the text gets
      
      if (6 <= length(levels(temp$Var1)) &
          length(levels(temp$Var1)) <= 10) {
        axis.size = 1
      }
      else if (length(levels(temp$Var1)) <= 6) {
        axis.size = 1.3
      }
      else if (20 <= length(levels(temp$Var1)) &
               length(levels(temp$Var1)) >= 10) {
        axis.size = .8
      } else {
        axis.size = .5
      }
      
      # margin text for the gene/variant names
      mtext(
        levels(temp$Var1), line = -49.5, at = mid.positions, adj = 1, las = 2, cex = 1
      )
      
      
      ### plot the axis
      # this is for the solid line
      
      axis(
        1, at = positions, las = 2, labels = rep(levels(temp$Subgroup),(max(positions) +
                                                                          1) / (length(
                                                                            levels(temp$Subgroup)
                                                                          ) + 1)), cex.axis = axis.size, hadj = 1
      )
      # this adds the tick marks below
      axis(
        1, at = c(0,remove.positions), las = 2, tck = -.075, lwd = 0, lwd.ticks =
          3, labels = NA, lty = 2
      )
      
      # tried to add segment lines, was a bit too hard to get them in the right place
      ### plot the lines under the axis
      #for(i in 1:length(seg.positions)){
      #  segments(seg.positions[[i]][1],-(par("usr")[4]-par("usr")[3])*(axis.size/5),seg.positions[[i]][2],-(par("usr")[4]-par("usr")[3])*(axis.size/5),  col = "#8B8B8B")
      #}
      ### add the text
      #text(x = mid.positions, y = par("usr")[3], labels = levels(temp$Var1), srt=90, cex = 1, adj = axis.size*4)
      
      
      
      #####
      #####
      ##### MULTI PLOTTER FROM HERE
      #####
      ####
      
      # 
      # this section plots all the individual variants on there own plot with
      # accompanying stats, these use the sink to keep that stat annotation
      # 
      
    } else {
      
      # turn on the sink
      LogSinker(
        data.type = data.type, outdir = paste0(output.dirname,"/Transcripts/", data.type, "/"), start = TRUE, cont = FALSE
      )
      
      #### multiple layout
      # this bit is for making the titles as above
      
      # was there a group.of.interest used?
      if (!any(group.of.interest == "all")) {
        # is yes then group.of.interest vs others (+/- the hgnc ID)
        if (is.character(hgnc.id[1,]) == "TRUE") {
          multi.title = paste0(
            hgnc.id," expression by ", data.type, " (",paste(group.of.interest, collapse  = ", "),") vs Others"
          )
        }
        else if (is.character(hgnc.id[1,]) == "FALSE") {
          multi.title = paste0(gene.of.interest," expression by ",data.type,  " (",paste(group.of.interest, collapse  =
                                                                              ", "),") vs Others")
        }
        
      } else {
        # no group.of.interest used, +/- hgnc.id
        if (is.character(hgnc.id[1,]) == "TRUE") {
          multi.title = paste0(hgnc.id," (", gene.of.interest, ") expression by ", data.type, " over subgroups")
        }
        else if (is.character(hgnc.id[1,]) == "FALSE") {
          multi.title = paste0(gene.of.interest," expression by ",data.type," over subgroups")
        }
        
      }
      
      
      
      # set amrgins
      par(mfrow = c(2,2), mar = c(5,6,5,6))
      
      ###plot as above with abline, this always plots the ENSG first, the rest are plotted below
      
      # first plot to be made is always the gene
      # xpd here as above for debugging
      par(xpd = FALSE)
      boxplot(
        temp$value[grep("^ENSG", temp$Var1)] ~ temp$Subgroup[grep("^ENSG", temp$Var1)], las = 1, col = sbgrp.col,
        ylab = "log2(FPKM+1)", main = gsub(data.type, gene.of.interest, title), ylim = scale
      )
      
      # draw a median line on
      abline(h = median(temp$value[grep("^ENSG", temp$Var1)]), lty = 2)
      
      # plot in margin
      par(xpd = TRUE)
      
      ### median text
      text(par("usr")[2], median(temp$value[grep("^ENSG", temp$Var1)]),"Median", cex = 1, adj = -.1)
      
      ### ANOVA
      # do this if there are more than two things to compare
      if (length(unique(temp$Subgroup[grep("^ENSG", temp$Var1)])) > 2) {
        # text for the log file
        cat("\nANOVA of ",  gsub(paste0(gene.of.interest, " \\(|) expression"), "", title),":\n\n")
        
        # do the anova and output the data for the sink
        aov(temp$value[grep("^ENSG", temp$Var1)] ~ temp$Subgroup[grep("^ENSG", temp$Var1)]) -> stat
        summary(stat) -> stat.summary
        print(stat.summary)
        
        # pull out the plotting values
        stat.summary[[1]]$"F value"[1] -> F
        stat.summary[[1]]$"Pr(>F)"[1] -> p.val
        
       # tidy them up 
        if (p.val < 0.001 & !is.na(p.val)) {
          p.text <- "p < 0.001"
        } else {
          p.text <- paste("p =" , round(p.val,3))
        }
        
        # write the stats on the chart
        f.text <- paste("F =",round(F,3))
        
        text(par("usr")[2],median(temp$value[grep("^ENSG", temp$Var1)]) - ((par("usr")[4] -
                                                                              par("usr")[3]) * .05), f.text, adj = -.1)
        text(par("usr")[2],median(temp$value[grep("^ENSG", temp$Var1)]) - ((par("usr")[4] -
                                                                              par("usr")[3]) * .1), p.text, adj = -.1)
        
        ### T TEST
      } else if (length(unique(temp$Subgroup[grep("^ENSG", temp$Var1)])) == 2) {
        # if there are onyl two things to compare
        # output for log
        cat("\nT Test of ",  gsub(paste0(gene.of.interest, " \\(|) expression"), "", title),":\n\n")
        
        # do the t test and print output
        t.test(temp$value[grep("^ENSG", temp$Var1)] ~ temp$Subgroup[grep("^ENSG", temp$Var1)]) -> stat
        
        print(summary(stat))
        
        # pull out the plotting values
        stat[[2]] -> df
        stat[[3]] -> p.val
        
        # tidy them up 
        if (p.val < 0.001 & !is.na(p.val)) {
          p.text <- "p < 0.001"
        } else {
          p.text <- paste("p =" , round(p.val,3))
        }
        
        f.text <- paste("DF =",signif(df,3))
        
        # write the stats on the chart
        text(par("usr")[2],median(temp$value[grep("^ENSG", temp$Var1)]) - ((par("usr")[4] -
                                                                              par("usr")[3]) * .05), f.text, adj = -.1)
        text(par("usr")[2],median(temp$value[grep("^ENSG", temp$Var1)]) - ((par("usr")[4] -
                                                                              par("usr")[3]) * .1), p.text, adj = -.1)
      }
      
      #### Get rid of the ENSG from above to plot the rest, get the names of the rest
      
      temp.noensg <- temp[grep("^ENSG",temp$Var1, invert = TRUE),]
      temp.names <- sort(levels(factor(temp.noensg$Var1)))
      
      cat("\nTrancript Names:\n\n")
      
      
      print(temp.names)
      #### loop over the names to plot as above with appropriate stats on them
      
      # plot all the rest of the variants, loop over them
      for (i in 1:length(temp.names)) {
        par(xpd = FALSE)
        # draw the plot
        boxplot(
          temp.noensg$value[grep(temp.names[i], temp.noensg$Var1)] ~  temp.noensg$Subgroup[grep(temp.names[i], temp.noensg$Var1)], col = sbgrp.col, main = gsub(data.type, temp.names[i], title),
          ylab = "log2(FPKM+1)", las = 1, ylim = scale
        )
        
        # median line
        abline(h = median(temp.noensg$value[grep(temp.names[i], temp.noensg$Var1)]), lty = 2)
        
        # allow to draw in margin
        par(xpd = TRUE)
        
        # median text
        text(
          par("usr")[2],median(temp.noensg$value[grep(temp.names[i], temp.noensg$Var1)]),"Median", cex = 1, adj = -.1
        )
        
        
        #### ANOVA
        
        if (length(unique(temp$Subgroup[grep("^ENSG", temp$Var1)])) > 2) {
          # more than two to compare? ANOVA
          # text for the sink
          cat("\nANOVA of ",  gsub(paste0(gene.of.interest, " \\(|) expression"), "", title),temp.names[i],":\n\n")
          
          # do an anova and print it for the sink
          aov(temp.noensg$value[grep(temp.names[i], temp.noensg$Var1)] ~ temp.noensg$Subgroup[grep(temp.names[i], temp.noensg$Var1)]) -> stat
          summary(stat) -> stat.summary
          
          print(stat.summary)
          
          # pull out and pretty the data
          
          stat.summary[[1]]$"F value"[1] -> F
          stat.summary[[1]]$"Pr(>F)"[1] -> p.val
          
          if (p.val < 0.001 & !is.na(p.val)) {
            p.text <- "p < 0.001"
          } else {
            p.text <- paste("p =" , round(p.val,3))
          }
          f.text <- paste("F =",round(F,3))
          
          # draw it on the chart
          text(par("usr")[2],median(temp.noensg$value[grep(temp.names[i], temp.noensg$Var1)]) -
                 ((par("usr")[4] - par("usr")[3]) * .05), f.text, adj = -.1)
          text(par("usr")[2],median(temp.noensg$value[grep(temp.names[i], temp.noensg$Var1)]) -
                 ((par("usr")[4] - par("usr")[3]) * .1), p.text, adj = -.1)
          
          #### T TEST
          
        } else if (length(unique(temp$Subgroup[grep("^ENSG", temp$Var1)])) == 2) {
          # for two variables do a mr T test, for the sink
          cat("\nT Test of ",  gsub(paste0(gene.of.interest, " \\(|) expression"), "", title),temp.names[i],":\n\n")
          
          # do the test, pull the data out and pretty up
          t.test(temp.noensg$value[grep(temp.names[i], temp.noensg$Var1)] ~ temp.noensg$Subgroup[grep(temp.names[i], temp.noensg$Var1)]) -> stat
          
          stat[[2]] -> df
          stat[[3]] -> p.val
          
          print(summary(stat))
          
          if (p.val < 0.001 & !is.na(p.val)) {
            p.text <- "p < 0.001"
          } else {
            p.text <- paste("p =" , round(p.val,3))
          }
          
          f.text <- paste("DF =",signif(df,3))
          
          # put the stats on the chart
          text(par("usr")[2],median(temp.noensg$value[grep(temp.names[i], temp.noensg$Var1)]) -
                 ((par("usr")[4] - par("usr")[3]) * .05), f.text, adj = -.1)
          text(par("usr")[2],median(temp.noensg$value[grep(temp.names[i], temp.noensg$Var1)]) -
                 ((par("usr")[4] - par("usr")[3]) * .1), p.text, adj = -.1)
        }
      }
      
      # turn the sink off
      LogSinker(start = FALSE)
      
    }
  }


plotBarData <- function(data, subgroup, output.dirname) {
  

  # 
  # as implied by the descriptive name this one plots bars
  # in particular this will draw the relative variant level for each subgroup
  # coloured by different things
  #
  # data = data to be plotted, e.g. isoform, cds
  # subgroup = slot [[6]] Subgroup
  # keep.index = slot [[4]] Keep Index Data
  # output.dirname = output directory used for the sink in multi portion of the plotter
  #

  
  # this identifies what sort of data we are dealing with and creates a
  # data.type for that
  
  if (any(grep("^ENST", rownames(data)))) {
    data.type = "Isoform"
  } else if (any(grep("^P", rownames(data)))) {
    data.type = "CDS"
  }else if (any(grep("^TSS", rownames(data)))) {
    data.type = "TSS"
  }
  
  
  # trun on the sink
  LogSinker(
    data.type = data.type, outdir = paste0(output.dirname,"/Transcripts/", data.type, "/"), start = TRUE, cont = TRUE
  )
  
  
  ### reverse factor to make the plot prettier
  factor.rev <- factor(subgroup, levels = rev(levels(subgroup)))
  ### get the subgroup means and remove any unwanted data ( basically the gene id)
  apply(data,1,function(x) {
    tapply(x, factor.rev, mean)
  }) -> means.subgroup
  grep("ENSG",colnames(means.subgroup)) -> remove.index
  # are there any to remove?
  if (length(remove.index) > 0) {
    means.subgroup <- means.subgroup[,-(remove.index)]
  }
  
  ###get the proportions and transpose the data
  t(prop.table(means.subgroup,1)) -> proportions.subgroup
  
  # for the sink give the percentages
  cat("\nTranscript Precentages:\n\n")
  
  print(proportions.subgroup)
  
  # One figure in row 1 and two figures in row 2
  # row 1 is 1/3 the height of row 2
  # column 2 is 1/4 the width of the column 1
  ### sort out the plotting area
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE),
         widths = c(1,1), heights = c(1,6))
  par(mar = c(0,0,3,0))
  # this is sorting out the positions for the plotting
  # positions for the x axis
  x.positions <- rep(c(0.2,0.5,0.8),6)
  # position for the numbers
  num.positions <- length(rownames(proportions.subgroup))
  # the number of rows
  ceiling(num.positions / 3) -> num.rows
  # y positions
  (1:num.rows) - 0.5 -> y.positions
  rep(y.positions, each = 3)[1:num.positions] -> y.positions
  
  ###colours
  # use the previous function to make colours
  assignColors(as.factor(rownames(proportions.subgroup))) -> colors
  
  ### plot the name of the isoform etc. in its colour
  try(plot(
    c(0,0.2,0.4,0.6,0.8,1),c(0,0,0,0,0,ceiling(num.positions / 3)), type = "n", axes = FALSE, ylab = "", xlab = "",  cex.main =
      2,   main = paste0(
        "Relative transcript levels by subgroup (coloured by: ", data.type,")"
      )
  ))
  text(
    x.positions[1:num.positions],y.positions, rownames(proportions.subgroup),
    col = colors, font = 2, cex = 1
  )
  ### plot the bars themselves
  par(mar = c(0,4,4,2))
  try(barplot(
    proportions.subgroup, las = 1,
    xlab = "Proportion of Transcript", col = colors,
    #legend = rownames(proportions.subgroup)
    axes = FALSE, horiz = TRUE, width = 0.2, ylab = ""
  ))
  axis(
    3, at = c(0,0.2,0.4,0.6,0.8,1),labels = c("0%","20%","40%","60%","80%","100%")
  )
  # turn off the sink
  LogSinker(start = FALSE)
}



plotTranscripts <-  function(gene.of.interest, color.by = "transcript", variants = TRUE) {
    
  # 
  # this function uses gViz to plot transcripts, can colour by transcript, tss
  # and cds also can add on variants to the bottom of the track the actual
  # variant table is plotted in the .Rmd
  # 
  # gene.of.interest = gene of interest
  # color.by = what to colour by...options are 'transcript', 'tss' and 'cds'
  # variants =  whether to plot the variants on too
  #
  #
  
  
  
  ### get the gtf together for the gene.of.interest
    # sort out the gene.of.interest/trailing .1 CAN THIS BE MOVED?
    gsub('\\.[0-9]+$', '', isoforms.attr$nearest_ref_id) -> ensembl.ids
  # get your gene id 
    as.character(isoforms.attr[which(ensembl.ids == gene.of.interest),"gene_id"]) -> gene.of.interest.id
    # split the gtf
    gtf[which(gtf$gene_id %in% gene.of.interest.id),] -> gtf.gene.of.interest
    # get rid of the ensembl ID
    grep("ENSG",gtf.gene.of.interest$oId) -> remove
    gtf.gene.of.interest[-remove,] -> gtf.gene.of.interest
    ###chromosome
    chr <- as.character(unique(seqnames(gtf.gene.of.interest)))
    ###genome
    # THIS IS HARDCODED BUT MAY NEED CHANGING IN THE FUTURE?
    gen <- "hg19"
    ### ideogram from above parameters with display settings
    itrack <- IdeogramTrack(genome = gen, chromosome = chr)
    displayPars(itrack) <- list(size = 0.5, littleTicks = TRUE)
    ### where it is on the genome plus pars
    gtrack <- GenomeAxisTrack()
    displayPars(gtrack) <- list(size = 1)
    
    #### sort out the varying tracks for tss/cds/isoform
    symbol.names <- paste(
      gsub('\\.[0-9]+$', '', gtf.gene.of.interest$nearest_ref),
      gtf.gene.of.interest$p_id,
      gtf.gene.of.interest$tss_id,
      gtf.gene.of.interest$gene_name,sep = "\n"
    )
    
    ### gr will be the only one plotted the others are used for colouring in the transcripts
    
    grtrack <-
      GeneRegionTrack(
        gtf.gene.of.interest, genome = gen,chromosome = chr, name = gene.of.interest,
        gene = gtf.gene.of.interest$gene_id,exon = gtf.gene.of.interest$exon_number,
        transcript = gtf.gene.of.interest$transcript_id, symbol =
          symbol.names
      )
    cdstrack <-
      GeneRegionTrack(
        gtf.gene.of.interest, genome = gen,chromosome = chr, name = gene.of.interest,
        gene = gtf.gene.of.interest$gene_id,exon = gtf.gene.of.interest$exon_number,
        transcript = gtf.gene.of.interest$p_id, symbol = symbol.names
      )
    tsstrack <-
      GeneRegionTrack(
        gtf.gene.of.interest, genome = gen,chromosome = chr, name = gene.of.interest,
        gene = gtf.gene.of.interest$gene_id,exon = gtf.gene.of.interest$exon_number,
        transcript = gtf.gene.of.interest$tss_id, symbol = symbol.names
      )
    
    ###colour them in according to which is being plotted
    
    if (color.by == "transcript") {
      fill.cols <- "#B3CDE3"
      number.groups <- nlevels(as.factor(group(grtrack)))
    }
    if (color.by == "tss") {
      assignColors(as.factor(group(tsstrack))) -> fill.cols
      number.groups <- nlevels(as.factor(group(tsstrack)))
    }
    if (color.by == "cds") {
      assignColors(as.factor(group(cdstrack))) -> fill.cols
      number.groups <- nlevels(as.factor(group(cdstrack)))
    }
    
    if (number.groups == 1) {
      font.size.group = 2
    }else{
      font.size.group = 1 / (log2(number.groups + 1) * 0.8)
    }
    
    ### sort out the display pars....still trying to work out how to put a title on the top!
    
    displayPars(grtrack) <-
      list(
        showId = TRUE, showExonId = FALSE, geneSymbol = TRUE,
        cex.group = font.size.group, fontcolor.title =
          "#808080",
        cex.title = 2,
        cex.ExonId = font.size.group,
        size = 1,
        fontcolor.item = "#808080",
        collapseTracks = FALSE,
        fill = as.character(fill.cols)
      )
    
    #### add in any variants
    if (variants == TRUE) {
      ###hardcoded for red and blue
      var.cols <-
        factor(c(rep("blue",nrow(
          dedup.variants.ns
        )),rep("red",nrow(
          dedup.variants.syn
        ))))
      
      ###variants generated outside of the function and fed in....do these need to be in the function?
      ### currently only coded for NS/SYN which will be a problem if later on people want to play with it.....
      ### rethink labelling of this
      
      temp.starts <-
        c(dedup.variants.ns$Start,dedup.variants.syn$Start)
      temp.ends <- c(dedup.variants.ns$End,dedup.variants.syn$End)
      no.features = nrow(dedup.variants.ns) + nrow(dedup.variants.syn)
      features <-
        factor(c(rep("NS",nrow(
          dedup.variants.ns
        )),rep("SYN",nrow(
          dedup.variants.syn
        ))), levels = c("NS", "SYN"))
      
      #### produce the output table for kable and also ids for later on IF they are present, switch from df to matrix to allow the same rowname
      
      ids <- c()
      out.df <- c()
      if (nrow(dedup.variants.ns) > 0) {
        id.ns <-
          paste0(
            dedup.variants.ns$sample,"\n",dedup.variants.ns$Chr,":",dedup.variants.ns$Start,"-",dedup.variants.ns$End,";\n",dedup.variants.ns$Ref,":",dedup.variants.ns$Alt
          )
        id.ns.df <-
          data.frame(
            dedup.variants.ns$sample,dedup.variants.ns$Chr,dedup.variants.ns$Start,dedup.variants.ns$End,paste(dedup.variants.ns$Ref,":",dedup.variants.ns$Alt)
          )
        colnames(id.ns.df) <-
          c("Sample","Chromosome", "Start", "End", "Variant")
        id.ns.df <- as.matrix(id.ns.df[order(id.ns.df$Start),])
        rownames(id.ns.df) <- rep("NS",nrow(id.ns.df))
        ids <- c(id.ns)
        out.df <- id.ns.df
      }
      if (nrow(dedup.variants.syn) > 0) {
        id.syn <-
          paste0(
            dedup.variants.syn$sample,"\n",dedup.variants.syn$Chr,":",dedup.variants.syn$Start,"-",dedup.variants.syn$End,";\n",dedup.variants.syn$Ref,":",dedup.variants.syn$Alt
          )
        id.syn.df <-
          data.frame(
            dedup.variants.syn$sample,dedup.variants.syn$Chr,dedup.variants.syn$Start,dedup.variants.syn$End,paste(dedup.variants.syn$Ref,":",dedup.variants.syn$Alt)
          )
        colnames(id.syn.df) <-
          c("Sample","Chromosome", "Start" , "End", "Variant")
        id.syn.df <- as.matrix(id.syn.df[order(id.syn.df$Start),])
        rownames(id.syn.df) <- rep("SYN",nrow(id.syn.df))
        ids <- c(ids,id.syn)
        out.df <- rbind(out.df, id.syn.df)
      }
      
      #### if no variants found do this
      
      if (nrow(dedup.variants.syn) + nrow(dedup.variants.ns) == 0) {
        plotTracks(
          list(itrack, gtrack,  grtrack), extend.left = 0.1,
          title.width = 2,cex.main = 5,
          background.panel = "white",background.title = "white",background.cex = 2, collapse = FALSE
        )
        
        
        return(NULL)
        
      } else {
        ### if they are make an annotation track, and plot that at the bottom
        
        antrack <-
          AnnotationTrack(
            start = temp.starts, end =  temp.starts,col.line = NULL,stacking = "full",
            group = features, feature = features, id = ids, genome = gen,chromosome = chr, name = "Variants"
          )
        
        displayPars(antrack) <-
          list(
            groupAnnotation = "group", min.width = 1.5,
            fill = as.character(var.cols), col = "#808080",
            cex = 1, col.line = "white", fontcolor.title =
              "#808080",
            cex.title = 2
          )
        
        plotTracks(
          list(itrack, gtrack,  grtrack, antrack), extend.left = 0.1,
          title.width = 2,cex.main = 5,
          background.panel = "white",background.title = "white",background.cex = 2, collapse = FALSE
        )
        
        return(out.df)
      }
      
      
    } else {
      ### plotting with no variants
      
      plotTracks(
        list(itrack, gtrack,  grtrack), extend.left = 0.1,
        title.width = 2,cex.main = 5,
        background.panel = "white",background.title = "white",background.cex = 2, collapse = FALSE
      )
      
      return(NULL)
      
    }
    
    
  }
