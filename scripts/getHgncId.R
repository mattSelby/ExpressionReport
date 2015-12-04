getHgncID <- function(gene.of.interest) {
  # this script is to create the gene symbol from the given ensembl ID and assigns to the envivronment
  # currently need to think about how the update from hg19 will change this
  #
  # gene.of.interest = ensembl ID for gene of interest
  #
  
  hgnc.id <-
    getBM(
      attributes = 'hgnc_symbol', filters = 'ensembl_gene_id', values = gene.of.interest, mart = ensembl
    )
  
  
  if (length(hgnc.id) > 0) {
    assign("hgnc.id", hgnc.id, .GlobalEnv)
    cat("HGNC ID Loaded\n")
    cat(paste("Gene:", as.character(hgnc.id)))
    return(hgnc.id)
    
  } else {
    cat("HGNC ID not Found\n")
  }
}