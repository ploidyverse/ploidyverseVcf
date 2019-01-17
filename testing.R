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

software(myvcf) <- c(Software = "polyRAD", Version = "1.1", Model = "IteratePopStruct",
                     Description = "Just some genotypes")
meta(header(myvcf))$ploidyverse
software(myvcf) <- rbind(software(myvcf),
                         DataFrame(row.names = "Imputation",
                                   Software = "Beagle",
                                   Version = "5000",
                                   Model = "spiffy polyploid",
                                   Description = "Imputation with Beagle"))

software(myvcf) <- meta(header(myvcf))$ploidyverse # takes vector, list, or table
software(myvcf)

software(myvcf) <- DataFrame(ID = c("GenotypeCalls", "Imputation"),
                             Software = c("polyRAD", "Beagle"),
                             Version = c("1.1", "5000"),
                             Model = c("IteratePopStruct", "Something"),
                             Description = c("Genotype calls from polyRAD",
                                             "Imputation with Beagle"))
software(myvcf)

software(myvcf) <- software(myvcf)

test <- DataFrame(A = c("a", "b", "c"), N = c(1,2,3))
test
row.names(test)
rownames(test) <- c('a', 'a', 'a')
test
