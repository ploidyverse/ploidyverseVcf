## Utility Functions to Help Format Data ##

# Take a 3D array of likelihood or probabilitiy values, in polyRAD format
# (allele copy number * sample * locus) and create a matrix formatted for
# VariantAnnotation (locus * sample, with a vector of values in each cell.)
array3D_to_matrixList <- function(arr){
  if(length(dim(arr)) != 3) stop("Need 3D array.")
  ngen <- dim(arr)[1]
  nsam <- dim(arr)[2]
  nal <- dim(arr)[3]
  outlist <- lapply(1:(nsam * nal), 
                    function(i) arr[((i - 1) * ngen + 1):(i * ngen)])
  outmat <- matrix(outlist, nrow = nal, ncol = nsam, 
                   dimnames = dimnames(arr)[3:2], byrow = TRUE)
  return(outmat)
}

# Reverse the above function, going from a matrix-list to a 3D array.
matrixList_to_array3D <- function(mat){
  nsam <- ncol(mat)
  nal <- nrow(mat)
  ngen <- unique(sapply(mat, length))
  if(length(ngen) > 1) 
    stop("Number of values not consistent across all samples and loci.")
  
  outarr <- array(unlist(t(mat)), dim = c(ngen, nsam, nal),
                  dimnames = list(as.character((1:ngen) - 1), colnames(mat),
                                  rownames(mat)))
  return(outarr)
}
