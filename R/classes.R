## Class definitions for ploidyverse VCF

setClass("ploidyverseVCF", contains = "CollapsedVCF")

# validity function for ploidyverse VCF before genotype calling
.valid.ploidyverseVCF.precall <- function(object){
  # call the VariantAnnotation validity function
  issues <- .valid.VCF(object)
  # check that allele depths exist
  if(!"AD" %in% names(geno(object))){
    issues <- c(issues, "No AD field for allele depth.")
  }
  # check for contig information
  # check for sample information
  return(issues)
}

# validity function for ploidyverse VCF after genotype calling
.valid.ploidyverseVCF.postcall <- function(object){
  issues <- .valid.ploidyverseVCF.precall(object)
  # check that genotype posterior probabilities exist
  if(!"GP" %in% names(geno(object))){
    issues <- c(issues, "No GP field for genotype probability")
  }
  return(issues)
}

# validity function for archival quality ploidyverse VCF
# check for tag sequences