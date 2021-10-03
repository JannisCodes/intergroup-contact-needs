BinaryCor <- function(x, type=c("pearson", "spearman"), ...){
  
  x <- as.matrix(x)
  C <- crossprod(x)
  R <- rcorr(x, type=type)$r
  p <- rcorr(x, type=type)$P
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste("N = ", diag(C), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[lower.tri(Rnew, diag = FALSE)] <- C[lower.tri(C, diag = FALSE)]
  Rnew <- as.data.frame(Rnew)
  
  ## remove last column and return the matrix (which is now a data frame)
  #Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew)
}