# Check that a VCF object is valid for various ploidyverse applications,
# and add appropriate header lines.
# These will look like:
# ##ploidyverseValidity=<ID=ploidyversePrecall,Valid=1,Description="File valid for calling genotypes with ploidyverse software.">
# ##ploidyverseValidity=<ID=ploidyversePostcall,Valid=1,Description="File contains genotype calls from ploidyverse software.">
# ##ploidyverseValidity=<ID=ploidyverseArchival,Valid=0,Description="File does not meet ploidyverse standards for data archiving.">
# where 1 indicates true and 0 indicates false.

setGeneric("markValidity", 
           function(object) standardGeneric("markValidity"))
setMethod("markValidity", "VCF", function(object){
  validout <- DataFrame(row.names = c("ploidyversePrecall", 
                                      "ploidyversePostcall",
                                      "ploidyverseArchival"),
                        Valid=rep(1L, 3))
  
  hdr <- header(object)
  
  # confirm allele depth present
  if(!"AD" %in% rownames(geno(hdr))){
    validout$Valid <- rep(0L, 3)
  } ## do any other checks to make sure it is formatted correctly?
  
  # confirm posterior probabilties exported
  if(!"GP" %in% rownames(geno(hdr))){
    validout$Valid[2] <- 0L
  }
  # confirm tag sequences provided if non-reference pipeline used
  if("NonRef" %in% seqnames(seqinfo(object)) && 
     !"Tag" %in% rownames(info(hdr))){
    validout$Valid[3] <- 0L
  }
  # confirm all contig info provided if reference pipeline used
  if(!"NonRef" %in% seqnames(seqinfo(object)) &&
     (any(is.na(seqlengths(seqinfo(object)))) || 
      any(is.na(genome(seqinfo(object)))))){
    validout$Valid[3] <- 0L
  }
  # confirm sample information provided
  if(!"META" %in% names(meta(hdr)) ||
     !"SAMPLE" %in% names(meta(hdr)) ||
     !"Species" %in% rownames(meta(hdr)$META) ||
     !"Ploidy" %in% rownames(meta(hdr)$META) ||
     !all(samples(hdr) %in% rownames(meta(hdr)$SAMPLE))){
    validout$Valid[3] <- 0L
  }
  
  # Add descriptions to the table
  desc <- rep("", 3)
  if(validout$Valid[1] == 1){
    desc[1] <- "File valid for calling genotypes with ploidyverse software."
  } else {
    desc[1] <- "File not valid for calling genotypes with ploidyverse software."
  }
  if(validout$Valid[2] == 1){
    desc[2] <- "File contains genotype calls from ploidyverse software."
  } else {
    desc[2] <- "File does not contain genotype calls from ploidyverse software."
  }
  if(validout$Valid[3] == 1){
    desc[3] <- "File meets ploidyverse standards for data archiving."
  } else {
    desc[3] <- "File does not meet ploidyverse standards for data archiving."
  }
  validout$Description <- desc
  
  # Add to VCF and return
  meta(header(object))$ploidyverseValidity <- validout
  return(object)
})

# Functions to return TRUE or FALSE to indicate whether a ploidyverse VCF
# is valid for various purposes.
# It may be wise to always run `markValidity` before running these functions,
# to ensure that nothing was modified that would have compromised validity.

.checkvalid <- function(object, type){
  # first confirm that ploidyverse validity is marked in the header
  if(!"ploidyverseValidity" %in% names(meta(header(object)))) return(FALSE)
  rname <- paste("ploidyverse", type, sep = "")
  if(!rname %in% 
     rownames(meta(header(object))$ploidyverseValidity)) return(FALSE)
  if(!"Valid" %in% 
     colnames(meta(header(object))$ploidyverseValidity)) return(FALSE)
  
  # then get the validity value
  v <- meta(header(object))$ploidyverseValidity[rname, "Valid"]
  return(v == 1 || v == "1")
}

setGeneric("validPloidyverseVCF_Precall",
           function(object) standardGeneric("validPloidyverseVCF_Precall"))
setMethod("validPloidyverseVCF_Precall", "VCF", function(object){
  return(.checkvalid(object, "Precall"))
})

setGeneric("validPloidyverseVCF_Postcall",
           function(object) standardGeneric("validPloidyverseVCF_Postcall"))
setMethod("validPloidyverseVCF_Postcall", "VCF", function(object){
  return(.checkvalid(object, "Postcall"))
})

setGeneric("validPloidyverseVCF_Archival",
           function(object) standardGeneric("validPloidyverseVCF_Archival"))
setMethod("validPloidyverseVCF_Archival", "VCF", function(object){
  return(.checkvalid(object, "Archival"))
})
