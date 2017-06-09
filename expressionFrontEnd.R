ExpressionReport <- function(gene.of.interest, output.dir, output.name = NULL, expression.preferences, expression.script) {
  
  Sys.time() -> TotTa
    require(knitr)
  require(xtable)
    require(rmarkdown)
    require(markdown)
  require(org.Hs.eg.db)
    source(expression.preferences)
  
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
  

 temp_gene_id  <- convertIDs(ids = gene.of.interest, fromKey = "ENSEMBL",toKey = "SYMBOL" , db = org.Hs.eg.db
             ,ifMultiple="useFirst")
  
    assign("gene.of.interest", gene.of.interest, envir = globalenv())
    output.dirname <- paste0(output.dir, temp_gene_id, "_",gene.of.interest,output.name,"/")
    assign("output.dirname", output.dirname, envir = globalenv())
    
    render(
      input = expression.script, 
      envir = globalenv(),
      output_file = paste0(temp_gene_id, "_",gene.of.interest,"_",  output.name, "_", "ExpressionReport.pdf"), 
      output_dir = output.dirname, 
      clean = TRUE)
    Sys.time() -> TotTb
    message(paste0("\nTotal Time: ",round(difftime(TotTb,TotTa, units = 'mins'),2),"min\n"))
  }


  