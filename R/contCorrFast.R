contCorrFast <- function(x,y,w, method=c("Pearson", "Spearman")) {
  if(!is.numeric(x)) { 
    x <- as.numeric(x)
  }
  if(!is.numeric(y)) {
    y <- as.numeric(y)
  }
    if(!is.numeric(w)) { 
    w <- as.numeric(w)
  }
  if(tolower(method[[1]])=="spearman") {
    x <- as.vector(wrankFast(x,w))
    y <- as.vector(wrankFast(y,w))
  }
  cont(x,y,w)
}
