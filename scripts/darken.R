darken <- function(color, factor=1.3){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  return(col)
}
