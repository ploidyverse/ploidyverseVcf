## Functions for converting genotype probabilities between multiallelic and
## quasi-biallelic format.

# Conversion matrix indicating how multiallelic genotypes correspond to 
# allele copy numbers.  Can do inverse if needing to convert from allele copy
# number probabilities to multiallelic genotype probabilities.
genoConvMat <- function(ploidy, n_alleles, inverse = FALSE){
  geno <- enumerateGenotypes(ploidy, n_alleles)
  outmat <- matrix(0L, nrow = (ploidy + 1) * n_alleles,
                   ncol = nrow(geno),
                   dimnames = list(paste(rep(0:(n_alleles - 1), each = ploidy + 1),
                                         rep(0:ploidy, times = n_alleles)),
                                   genotypeStrings(ploidy, n_alleles, sep = "")))
  
  for(g in 1:nrow(geno)){
    for(a in 1:n_alleles){
      thisrow <- (a - 1) * (ploidy + 1) + geno[g, a] + 1
      outmat[thisrow, g] <- 1L
    }
  }
  
  if(inverse){
    invmat <- MASS::ginv(outmat)
    rownames(invmat) <- colnames(outmat)
    colnames(invmat) <- rownames(outmat)
    return(invmat)
  } else {
    return(outmat)
  }
}

# Take a 3D array of allele copy number probabilities and convert to 
# multiallelic genotype probabilities.  3D array is formatted as in polyRAD.
acn_to_geno <- function(probarray, alleles2loc){
  ploidyp1 <- dim(probarray)[1]
  ploidy <- ploidyp1 - 1
  nind <- dim(probarray)[2]
  
  # get all numbers of alleles per locus
  al_n <- sort(unique(table(alleles2loc)))
  # get all conversion matrices
  conv_list <- lapply(al_n, function(n) genoConvMat(ploidy, n, inverse = TRUE))
  
  # matrix-list to output, loci x samples
  outmat <- matrix(list(), nrow = max(alleles2loc), ncol = dim(probarray)[2])
  
  # loop through loci and fill in the matrix
  for(L in unique(alleles2loc)){
    nal <- dim(thisarr)[3]
    thisarr <- probarray[,, alleles2loc == L]
    thisAlMat <- matrix(aperm(thisarr, c(1, 3, 2)),
                        nrow = ploidyp1 * nal,
                        ncol = nind)
    thisGenMat <- conv_list[[match(nal, al_n)]] %*% thisAlMat
    outmat[L,] <- lapply(1:nind, function(i) unname(thisGenMat[,i]))
  }
  
  return(outmat)
}
