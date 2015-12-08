


#### read in fusion table and sort names out
#mb_fusions <-  read.delim("/home/data/pbt/RNASeq/report_data/MB_20_April_2015/fusions/all_MB_results.txt" )
#short.sample.name <- gsub("Sample_","",mb_fusions$sample.name)
#short.sample.name <- gsub("_","",short.sample.name)
#### create just NMB names and filter for those on exprs data
#mb_fusions <- cbind(mb.fusions,short.sample.name)
#mb.fusions.short <-  mb.fusions[mb.fusions$short.sample.name  %in% annot$Index,]


### function to make fusion table
fusionMaker <- function(mb.fusions.short = mb.fusions.short, hgnc.id = hgnc.id, cutoff = cutoff, SOAPdir = "/home/dan/SOAP/output", output.dirname = output.dirname) {
  
  #
  # This function is to pull out any fusions that have been found for the sample and dumps all the .svg 
  # SOAP images into the output file
  #
  # mb.fusions.short = the fusion list generated from the load data function
  # hgnc.id = the gene id created earlier
  # cutoff = the minmum number of spanning reads created by SOAP for it to be classed as a fusion
  # SOAPdir = this is where the fusions live
  # THIS NEEDS TO BE MOVED TO THE CENTRAL DATA STORE
  # output.dirname = location for the log file
  #
  #

  # start the sink and write the cutoff
  
    LogSinker(
      data.type = "Fusions", outdir = paste0(output.dirname,"/Fusions/"), start = TRUE
    )
    
    cat(paste0("\nCutoff: ", cutoff, "\n\n"))
    
    #### search up and down gene for gene.of.interest and filter on cutoff
    # checking both sides of the fusion genes to get all possible fusions
    # using ^ and $ to get exact matches of the gene names
    #
    
    gene.of.interest.fusion <-
      mb.fusions.short[with(
        mb.fusions.short, grepl(paste0("^",hgnc.id,"$"), mb.fusions.short$up_gene) |
          grepl(paste0("^",hgnc.id,"$"), mb.fusions.short$dw_gene)
      ),]
    
    # this then checks that there are higher reads in BOTH the spanning and junction reads compared to the cutoff
    
    gene.of.interest.fusion <-
      gene.of.interest.fusion[gene.of.interest.fusion$Span_reads_num >= cutoff &
                   gene.of.interest.fusion$Junc_reads_num >= cutoff,]
    
    # print the table for the log and only proceed if more than one row
    print(gene.of.interest.fusion)
    if (nrow(gene.of.interest.fusion) > 0) {
      ### data frame for the table
      # looks complicated but is only prettying up the text so it looks sensible when you put it into the kable
      gene.of.interest.fusion.small <-
        data.frame(
          paste(gene.of.interest.fusion$up_gene, ":", gene.of.interest.fusion$dw_gene),
          paste0(
            gene.of.interest.fusion$up_chr," (",gene.of.interest.fusion$up_Genome_pos, ") : ", gene.of.interest.fusion$dw_chr, " (",gene.of.interest.fusion$dw_Genome_pos, ")"
          ),
          paste0(
            gene.of.interest.fusion$Span_reads_num, "/", gene.of.interest.fusion$Junc_reads_num
          ),
          paste0(gene.of.interest.fusion$up_loc, "/", gene.of.interest.fusion$dw_loc),
          paste0(gene.of.interest.fusion$up_strand, "/", gene.of.interest.fusion$dw_strand),
          paste0(
            gene.of.interest.fusion$Fusion_Type, " (", gene.of.interest.fusion$down_fusion_part_frame.shift_or_not, ")"
          )
        )
      # put into a matrix rather than data frame so we can have the same row names
      gene.of.interest.fusion.small <- as.matrix(gene.of.interest.fusion.small)
      colnames(gene.of.interest.fusion.small) <-
        c(
          "Fusion", "Chromosome","Span/Junction Reads","Up/Down Locus" ,"Up/Down Strand", "Type"
        )
      rownames(gene.of.interest.fusion.small) <- gene.of.interest.fusion$short.sample.name
      
      # order to make it pretty
      gene.of.interest.fusion.small <- gene.of.interest.fusion.small[order(rownames(gene.of.interest.fusion.small)),]
      
      
      ### output list
      ### [[1]] table
      ### [[2]] search param for next function, it has the search terms to find the svg files
      ### [[3]] directory list for searching, directory locations for the next script
      fusion.output.list <- list()
      fusion.output.list[[1]] <- gene.of.interest.fusion.small
      fusion.output.list[[2]] <- unique(paste0("/", gene.of.interest.fusion$up_gene, "_", gene.of.interest.fusion$dw_gene,"/"))
      fusion.output.list[[3]] <- dir(SOAPdir)[grep(paste(unique(gene.of.interest.fusion$sample_name),collapse = "|"), dir(SOAPdir))]
      
      # turn off the sink and return the list
      LogSinker(start = FALSE)
      
      return(fusion.output.list)
    }
  }


unTarGetSvg <- function(sample.name = "Sample_NMB166.tar.gz", SOAPdir = "/home/dan/SOAP/output",outdir = "/home/shelby/report/fusion.output", search.params = search.params) {
    
  # this script only runs if there are fusions which were identified from the search.params
  #
  # function is to untar and pull the svg of the fusions identified and transfer them to the output directory
  #
  # sample.name =  this comes from the fusion_maker script it contains all the smaples which were identified as having fusions
  # looped over in the .Rmd
  # SOAPdir = this is where the fusions live
  # output.dirname = location for the log file
  # search.params =  the search criteria for each of the samples which were dientified as interesting
  #
  
  # this will onyl run if there were fusions found
  ### are there fusions?
    if (any(!search.params == "/_/"))  {
      
      # this is to just re set the wd at the end of the function
      current.wd <- getwd()
      # list all files in the SOAP dir, find the corresponding names and untar, the long string from the search params only pulls out unique things
      setwd(SOAPdir)
      temp.files <-
        untar(tarfile = sample.name, list = TRUE)[grep(paste0("(?=.*",search.params,")(?=.*.svg)",collapse = "|"), perl = TRUE, untar(tarfile = sample.name, list = TRUE))]
      
      # untar the above only search the names
      untar(tarfile = sample.name, files = temp.files, exdir = paste0(outdir,"/",sample.name))
      
      #### move the files to clean everything up and delete
      file.rename(from = paste0(outdir,"/",sample.name,"/",temp.files), to = paste0(
          outdir,"/",sample.name,"/",gsub("\\S+/","", temp.files,perl = TRUE))
      )
      
      # delete the long files after moving them up levelst
      unlink(paste0(outdir,"/",sample.name,"/home"), recursive = TRUE)
      
      #out.dir<-paste0(outdir,"/",sample.name)
      #setwd(paste0(outdir,"/",sample.name))
      #name.list<-list.files(pattern = "*.svg")
      #for ( i in 1:length(name.list)){
      #eval(parse( text = paste("system(paste(c('convert', name.list[i], gsub('.svg', '.png', name.list[i])), collapse = ' '))")))
      
      #}
      # return to normal wd
      setwd(current.wd)
      
    }
    
  }
