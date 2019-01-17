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
  if(!all(sam %in% rownames(value))){
    warning("Some samples in VCF not found in row names of table")
    cat(sam[!sam %in% rownames(value)], sep = "\n")
  }
  if("factor" %in% unlist(lapply(value, class))){
    warning("Factor columns found; should they be character?")
  }
  if(!all(c("Species", "Ploidy") %in% colnames(value))){
    warning("Both Species and Ploidy columns will need to be present to meet ploidyverse archival specifications.")
  }
  if("Ploidy" %in% colnames(value)){
    value$Ploidy <- trimws(as.character(value$Ploidy))
    if(any(!grepl("^[[:digit:]]x(\\+[[:digit:]]x)*$", value$Ploidy))){
      stop("Ploidy not formatted correctly.")
    }
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
  if(any(!colnames(value) %in% rownames(metatab))){
    extcol <- colnames(value)[!colnames(value) %in% rownames(metatab)]
    metatab <- rbind(metatab, DataFrame(row.names = extcol,
                                        Type = rep("String", length(extcol)),
                                        Number = rep(".", length(extcol)),
                                        Description = rep("Insert description here", length(extcol))))
    warning("Non-standard column names found; please add descriptions to meta(header(object))$META")
  }
  # subset to only columns used
  metatab <- metatab[colnames(value),]
  
  # modify VCF and return
  meta(header(object))$SAMPLE <- DataFrame(value)
  meta(header(object))$META <- metatab
  return(object)
})

# specify software version
setGeneric("software<-",
           function(object, value) standardGeneric("software<-"))
setReplaceMethod("software", "VCF", function(object, value){
  if(is.null(names(value)) ||
     !all(c("Software", "Version", "Model", "Description") %in% names(value))){
    stop("Replacement value must be named list or vector with names Software, Version, Model, and Description.")
  }
  if(is(value, "DataFrame")){
    softtable <- value
  } else {
    softtable <- do.call(DataFrame, args = as.list(value))
  }
  if("ID" %in% colnames(softtable)){
    rownames(softtable) <- softtable$ID
    softtable <- softtable[,-match("ID", colnames(softtable))]
  } else {
    if(nrow(softtable) > 1 && (is.null(rownames(softtable)) || 
                               length(unique(rownames(softtable))) < 
                                 nrow(softtable))){
      stop("Multiple entries; please provide unique ID field.")
    }
    if(nrow(softtable) == 1 && is.null(rownames(softtable))){
      rownames(softtable) <- "GenotypeCalls"
    }
  }
  
  meta(header(object))$ploidyverse <- softtable ## maybe change name of this table
  return(object)
})
setGeneric("software", function(object) standardGeneric("software"))
setMethod("software", "VCF", function(object){
  if(is.null(meta(header(object))$ploidyverse)){
    return(DataFrame(Software = character(0),
                     Version = character(0),
                     Model = character(0),
                     Description = character(0)))
  } else {
    return(meta(header(object))$ploidyverse)
  }
})
