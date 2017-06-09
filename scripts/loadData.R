loadData <- function(db, load.variant.df = TRUE) {
  # All packages should be loaded in here.
  # Check to see all are installed this will screw up the knitr if not
  #
  # db = location of the folders with the data in
  # load.variant.df = whether to create the variant/grange objects(FALSE) or load in premade ones (TRUE)
  #
  
  library(biomaRt)
  library(reshape2)
  library(minfi)
  library(rtracklayer)
  library(knitr)
  library(Gviz)
  library(GenomicRanges)
  library(survival)
  library(scales)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  #source("~/report/script/continue.script.R")
  
  # for each of the data sets search in the global environment, if it's not there load it in
  # as a note this is set to automatically load the data into the working environment
  
  # Ensembl load for gene names
  #if (length(grep("^ensembl$", ls(envir = .GlobalEnv))) == 1) {
  #  cat("Ensembl Pre Loaded\n")
  #} else {
  #  #ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
  #  
  #  # hardcoded to hg19, does this need changing somehow?
  #  ensembl=useMart(host='feb2014.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
  #  assign("ensembl", ensembl, .GlobalEnv)
  #  cat("Ensembl Loaded\n")
  #  
  #}
  
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
  

  
  
  # Annotation files for expression
  if (length(grep("^annot$", ls(envir = .GlobalEnv))) == 1) {
    cat("Annot Pre Loaded\n")
  } else {
    annot <- read.delim(paste0(db,"/samples_data/sample_data.txt"))
    assign("annot", annot, .GlobalEnv)
    cat("Annot Loaded\n")
    
  }
  
  # VSD for RNA seq
  if (length(grep("^vsd$", ls(envir = .GlobalEnv))) == 1) {
    cat("VSD Pre Loaded\n")
  } else {
    vsd <-
      read.delim(paste0(db,"/expression/primary/genes/vsd/vsd.mb.txt"))
    
    colnames(vsd) <- gsub("\\..*", "", colnames(vsd))
    assign("vsd", vsd, .GlobalEnv)
    cat("VSD Loaded\n")
    
  }
  
  
  # The three data sets for CDS (counts, fpkm and attr)
  if (length(grep("^cds.attr$", ls(envir = .GlobalEnv))) + length(grep("^cds.counts$", ls(envir =
                                                                                          .GlobalEnv))) + length(grep("^cds.fpkm$", ls(envir = .GlobalEnv))) == 3) {
    cat("CDS Pre Loaded\n")
  } else {
    cds.fpkm <-
      read.delim(paste0(db,"/expression/primary/CDS/cds.fpkm_table"), row.names = 1)
    cds.counts <-
      read.delim(paste0(db,"/expression/primary/CDS/cds.count_table"), row.names = 1)
    cds.attr <-
      read.delim(paste0(db,"/expression/primary/CDS/cds.attr_table"))
    assign("cds.fpkm", cds.fpkm, .GlobalEnv)
    assign("cds.counts", cds.counts, .GlobalEnv)
    assign("cds.attr", cds.attr, .GlobalEnv)
    cat("CDS Loaded\n")
  }
  
  # The three data sets for isoforms (counts, fpkm and attr)
  if (length(grep("^isoforms.attr$", ls(envir = .GlobalEnv))) + length(grep("^isoforms.counts$", ls(envir =
                                                                                                    .GlobalEnv))) + length(grep("^isoforms.fpkm$", ls(envir = .GlobalEnv))) ==
      3) {
    cat("Isoforms Pre Loaded\n")
  } else {
    isoforms.fpkm <-
      read.delim(paste0(db,"/expression/primary/isoforms/isoforms.fpkm_table"), row.names = 1)
    isoforms.counts <-
      read.delim(paste0(db,"/expression/primary/isoforms/isoforms.count_table"), row.names = 1)
    isoforms.attr <-
      read.delim(paste0(db,"/expression/primary/isoforms/isoforms.attr_table"))
    assign("isoforms.fpkm", isoforms.fpkm, .GlobalEnv)
    assign("isoforms.counts", isoforms.counts, .GlobalEnv)
    assign("isoforms.attr", isoforms.attr, .GlobalEnv)
    cat("Isoforms Loaded\n")
  }
  
  
  # The three data sets for tss (counts, fpkm and attr)
  if (length(grep("^tss.attr$", ls(envir = .GlobalEnv))) + length(grep("^tss.counts$", ls(envir =
                                                                                          .GlobalEnv))) + length(grep("^tss.fpkm$", ls(envir = .GlobalEnv))) == 3) {
    cat("TSS Pre Loaded\n")
  } else {
    tss.fpkm <-
      read.delim(paste0(db,"/expression/primary/TSS/tss_groups.fpkm_table"), row.names = 1)
    tss.counts <-
      read.delim(paste0(db,"/expression/primary/TSS/tss_groups.count_table"), row.names = 1)
    tss.attr <-
      read.delim(paste0(db,"/expression/primary/TSS/tss_groups.attr_table"))
    assign("tss.fpkm", tss.fpkm, .GlobalEnv)
    assign("tss.counts", tss.counts, .GlobalEnv)
    assign("tss.attr", tss.attr, .GlobalEnv)
    cat("TSS Loaded\n")
  }
  
  
  # The three data sets for genes (counts, fpkm and attr)
  if (length(grep("^genes.attr$", ls(envir = .GlobalEnv))) + length(grep("^genes.counts$", ls(envir =
                                                                                              .GlobalEnv))) + length(grep("^genes.fpkm$", ls(envir = .GlobalEnv))) == 3) {
    cat("Genes Pre Loaded\n")
  } else {
    genes.fpkm <-
      read.delim(paste0(db,"/expression/primary/genes/genes.fpkm_table"), row.names = 1)
    genes.counts <-
      read.delim(paste0(db,"/expression/primary/genes/genes.count_table"), row.names = 1)
    genes.attr <-
      read.delim(paste0(db,"/expression/primary/genes/genes.attr_table"))
    assign("genes.fpkm", genes.fpkm, .GlobalEnv)
    assign("genes.counts", genes.counts, .GlobalEnv)
    assign("genes.attr", genes.attr, .GlobalEnv)
    cat("Genes Loaded\n")
  }
  
  #load in the gencode gtf
  if (length(grep("^gtf$", ls(envir = .GlobalEnv))) == 1) {
    cat("GTF Pre Loaded\n")
  } else {
    gtf <-
      import("/home/dan/cuffdiff_results/Gencode_v17_Matt.combined.gtf", format = "gtf")
    assign("gtf", gtf, .GlobalEnv)
    cat("GTF Loaded\n")
  }
  
  ###THIS NEEDS FIXING
  #if(length(grep("^meth_primary_detP$", ls(envir=.GlobalEnv)))+length(grep("^meth_primary$", ls(envir=.GlobalEnv)))==2
  # ) { cat("Meth Pre Loaded\n")
  #} else {
  #  meth_primary_detP <- readRDS(paste0(db,"/methylation/primary/MB_450k_Primary_DetPValues.rds"))
  # meth_primary <- readRDS(paste0(db,"/methylation/primary/MB_450k_Primary.rds"))
  # assign("meth_primary_detP", meth_primary_detP, .GlobalEnv)
  # assign("meth_primary", meth_primary, .GlobalEnv)
  # cat("Meth Loaded\n")
  #}
  
  # this is a bit more complicated, on updating the variant.df must be created but this can then be save and loaded from here
  # won't work overly well with the current loading system...
  # also create a grange object
  # to create these the script is hashed below
  # NEED TO FIX THIS
  if (length(grep("^gencode.gtf$", ls(envir = .GlobalEnv))) + length(grep("^nmb.all.vars.filt$", ls(envir =
                                                                                                    .GlobalEnv))) + length(grep("^variants.df$", ls(envir = .GlobalEnv))) ==
      3)
  {
    cat("Variants Pre Loaded\n")
  } else {
    ### sort out the names of the variants here
    load(paste0(db,"/mutations/nmb.all.vars.filt.Rdata"))
    nmb.all.vars.filt$sample <-
      gsub("Sample_","",nmb.all.vars.filt$sample)
    nmb.all.vars.filt$sample <-
      gsub("_","",nmb.all.vars.filt$sample)
    assign("nmb.all.vars.filt", nmb.all.vars.filt, .GlobalEnv)
    gencode.gtf <-
      import("/home/dan/gencode.v17.annotation.exon.cds.gtf", format = "gtf")
    assign("gencode.gtf", gencode.gtf, .GlobalEnv)
    
    if (load.variant.df == TRUE) {
      load(file = "/home/shelby/report/variants.df.RData")
      assign("variants.df", variants.df, .GlobalEnv)
      load(file = "/home/shelby/report/grange.all.vars.RData")
      assign("grange.all.vars", grange.all.vars, .GlobalEnv)
    } else {
      ### option to remake this variant table and grange
      variants.df <-
        data.frame(
          sampleNames = as.character(nmb.all.vars.filt$sample), seqnames = as.factor(as.character(nmb.all.vars.filt$Chr)),
          starts = as.numeric(as.character(nmb.all.vars.filt$Start)) -
            1,
          ends = as.numeric(as.character(nmb.all.vars.filt$End)),
          names = c(rep(".", nrow(
            nmb.all.vars.filt
          ))),
          scores = c(rep(".", nrow(
            nmb.all.vars.filt
          ))),
          strands = c(rep("+", nrow(
            nmb.all.vars.filt
          )))
        )
      
      #save(variants.df, file = "variants.df.RData")
      
      assign("variants.df", variants.df, .GlobalEnv)
      
      grange.all.vars <- GRanges(
        seqnames = variants.df$seqnames,
        ranges = IRanges(as.numeric(as.character(
          variants.df$starts
        )),
        as.numeric(as.character(
          variants.df$ends
        ))),
        strand = variants.df$strands
      )
      #save(grange.all.vars, file = "grange.all.vars.RData")
      
      assign("grange.all.vars", grange.all.vars, .GlobalEnv)
      
    }
    
    cat("Variants Loaded\n")
  }
  
  
  #load in the allele specific expression data
  if (length(grep("^ase.data$", ls(envir = .GlobalEnv))) == 1) {
    cat("ASE Pre Loaded\n")
  } else {
    load(paste0(db,"/expression/primary/ASE/aseData.Rdata"))
    assign("ase.data", ase.data, .GlobalEnv)
    cat("ASE Loaded\n")
    
  }
  
  
  #load in the fusion data
  if (length(grep("^mb.fusions.short$", ls(envir = .GlobalEnv))) == 1) {
    cat("Fusions Pre Loaded\n")
  } else {
    
    #### read in fusion table and sort names out
    mb.fusions <-  read.delim(paste0(db,"/fusions/all_MB_results.txt"))
    short.sample.name <- gsub("Sample_","",mb.fusions$sample_name)
    short.sample.name <- gsub("_","",short.sample.name)
    ### create just NMB names and filter for those on exprs data
    mb.fusions <- cbind(mb.fusions,short.sample.name)
    mb.fusions.short <-  mb.fusions[mb.fusions$short.sample.name  %in% annot$Index,]
    
    assign("mb.fusions.short", mb.fusions.short, .GlobalEnv)
    cat("Fusions Loaded\n")
    
  }
  
  #load in the cufflinks data
  if (length(grep("^cuff.samp$", ls(envir = .GlobalEnv))) == 1) {
    cat("Cuff Pre Loaded\n")
  } else {
    cuff.samp <-  read.delim(paste0(db,"/expression/primary/cufflinks.samples.table"))
    assign("cuff.samp", cuff.samp, .GlobalEnv)
    cat("Cuff Loaded\n")
  }

  
  
}
