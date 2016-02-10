ExpressionReport <- function(gene.of.interest, output.dir, output.name = NULL, expression.preferences, expression.script) {
  
  Sys.time() -> TotTa
    require(knitr)
    require(rmarkdown)
    require(markdown)
    source(expression.preferences)
  
  
    assign("gene.of.interest", gene.of.interest, envir = globalenv())
    output.dirname <- paste0(output.dir,gene.of.interest,output.name,"/")
    assign("output.dirname", output.dirname, envir = globalenv())
    
    render(
      input = expression.script, envir = globalenv(),output_file = paste0(gene.of.interest,"_", output.name, "_", "ExpressionReport.pdf"), 
      output_dir = output.dirname, clean = FALSE)
    Sys.time() -> TotTb
    message(paste0("\nTotal Time: ",round(difftime(TotTb,TotTa, units = 'mins'),2),"min\n"))
  }


  