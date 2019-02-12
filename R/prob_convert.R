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
