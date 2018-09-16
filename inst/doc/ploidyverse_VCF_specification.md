# Variant Call Format specification for the ploidyverse

The Variant Call Format, or VCF, is the standard in bioinformatics for storing
information about DNA sequence variation.  One VCF file can contain genotypes
across many loci and individuals.  Because VCF is highly flexible, widely used,
and self-documenting, we adopt it here as the primary format for transferring
data among ploidyverse packages, and between the ploidyverse and other
software.

VCF is stored as tab-delimited text, with a corresponding binary format called
BCF.  Additionally, the 
[VariantAnnotation](https://doi.org/doi:10.18129/B9.bioc.VariantAnnotation) 
package from the Bioconductor 
collection of R packages implements flexible S4 classes for storing any or all
data from a VCF file within an R object, and has functions for reading and 
writing VCF files.  In the ploidyverse, we extend these S4 classes, adding
requirements for data that are mandatory for ploidyverse packages, and
generally setting standards for storing genotype data from polyploid organisms.

## Guiding principles

VCF data within the ploidyverse should:

* Be fully compatible with the 
[current VCF specification](https://github.com/samtools/hts-specs).
* Be self-documenting.  Someone 50 years from now should be able to open up
the file and understand what is going on and how to use the data.
* Be adapted for use with non-model organisms.  There may be no reference
genome, or a reference genome for a different species from the one being
studied.
* Provide allelic read depths and posterior probability estimates for all
genotypes.  Allele dosage will be uncertain without very high depth sequencing,
and that uncertainty must be quantified in the file.
