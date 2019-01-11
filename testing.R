library(VariantAnnotation)

myvcf <- 
  readVcf(system.file("extdata", "Msi01genes.vcf", package = "polyRAD"),
          param = ScanVcfParam(geno = "AD"))
meta(header(myvcf))
mysam <- samples(header(myvcf))
mysam
"PMS-430" %in% mysam

sampleinfo(myvcf) <- data.frame(row.names = mysam,
                                Species = rep("Miscanthus sinensis", length(mysam)),
                                Ploidy = rep("2x", length(mysam)),
                                stringsAsFactors = FALSE)
meta(header(myvcf))$META
sampleinfo(myvcf)
