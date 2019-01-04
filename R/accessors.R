# Accessors for some ploidyverse-specific information

setGeneric("sampleinfo",
           function(object) standardGeneric("sampleinfo"))
setMethod("sampleinfo", "VCF", function(object){
  return(meta(header(object))$SAMPLE)
})

setGeneric("sampleinfo<-",
           function(object, value) standardGeneric("sampleinfo<-"))
setReplaceMethod("sampleinfo", "VCF", function(object, value){
  sam <- samples(header(object))
  if(all(!rownames(value) %in% sam)){
    stop("Need row names that match sample names from VCF")
  }
  if(!all(rownames(value) %in% sam)){
    warning("Some row names not found in sample names from VCF")
    cat(rownames(value)[!rownames(value) %in% sam], sep = "\n")
  }
  if("factor" %in% unlist(lapply(value, class))){
    warning("Factor columns found; should they be character?")
  }
  
  # set up table of column information
  metatab <- DataFrame(row.names = c("Species", "Ploidy", "CollectionLocation",
                                     "MaterialProvider", "PopulationDesign",
                                     "BioSample", "DOI"),
                       Type = rep("String", 7), Number = rep(".", 7),
                       Description = c("Species name",
                                       "Ploidy with respect to reference genome",
                                       "Collection location of wild samples",
                                       "Person or organization curating biological materials",
                                       "Identity within artificial population",
                                       "Accession number in NCBI BioSample database",
                                       "Digital Object Identifier"))
  # check for any non-standard columns, print message about how to add description
  # subset to only columns used
  
  # modify VCF and return
  meta(header(object))$SAMPLE <- DataFrame(value)
  meta(header(object))$META <- metatab
  return(object)
})
