#### ask dan what each of these mean


filteringAnno <- function(x, depth = 50, perc = 0.3, filters.include = "all"){
  
    
    if(any(filters.include==(paste(filters.include,collapse = "|")))){
      filters.include <- c("X1000g.filt","snp.filt","cg46.filt","esp6500si_all.filt","depth.filt","perc.filt","NS.filt","poly.a.filt","non.exonic.filt","del.filt","filt.data.frame")
    }
    
  
  
  i=1
  filter.list <- list()
  if("X1000g.filt"%in%filters.include){
    x$X1000g2012apr_eur=="." -> X1000g.filt
    filter.list[[i]] <- X1000g.filt
    i+1 ->i
  }
  if("snp.filt"%in%filters.include){
    x$snp138=="." -> snp.filt
    filter.list[[i]] <- snp.filt
    i+1 ->i
  }
  if("cg46.filt"%in%filters.include){
    x$cg46=="." -> cg46.filt
    filter.list[[i]] <- cg46.filt
    i+1 ->i
  }
  if("esp6500si_all.filt"%in%filters.include){
    x$esp6500si_all=="." -> esp6500si_all.filt
    filter.list[[i]] <- esp6500si_all.filt
    i+1 ->i
  }
  if("depth.filt"%in%filters.include){
    x$total.depth > depth -> depth.filt
    filter.list[[i]] <- depth.filt
    i+1 ->i
  }
  if("perc.filt"%in%filters.include){
    x$perc.alt > perc -> perc.filt
    filter.list[[i]] <- perc.filt
    i+1 ->i
  }
  if("NS.filt"%in%filters.include){
    x$ExonicFunc.refGene=="nonsynonymous SNV" -> NS.filt
    filter.list[[i]] <- NS.filt
    i+1 ->i
  }
  if("poly.a.filt"%in%filters.include){
    (x$Alt%in%c("A","AA","AAA","AAAA","AAAAA","AAAAAA","AAAAAAA") & x$Ref == "-")|(x$Alt%in%c("T","TT","TTT","TTTT","TTTTT","TTTTTT","TTTTTTT") & x$Ref == "-") -> poly.a.filt
    filter.list[[i]] <- !poly.a.filt
    i+1 ->i
  }
  if("non.exonic.filt"%in%filters.include){
    x$Func.wgEncodeGencodeBasicV19 !="exonic"-> non.exonic.filt
    filter.list[[i]] <-  non.exonic.filt
    i+1 ->i
  }
  if("del.filt"%in%filters.include){
    x$LJB23_RadialSVM_pred == "D" -> del.filt
    filter.list[[i]] <- del.filt
    i+1 ->i
  }
  
  do.call(cbind, filter.list) -> filt.data.frame
  filt.data.frame[which(is.na(filt.data.frame))] <- FALSE
  apply(filt.data.frame,1,all) -> keep.index
  return(x[keep.index,])
}