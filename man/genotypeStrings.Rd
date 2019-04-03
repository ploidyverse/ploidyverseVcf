\name{genotypeStrings}
\alias{genotypeStrings}
\title{
Enumerate Strings for Representing All Genotypes in GT Field
}
\description{
For any given ploidy and number of alleles, the 
\href{https://github.com/samtools/hts-specs}{VCF specification} indicates the
order in which genotypes should be listed in the \code{GP}, \code{GL}, and
\code{PL} fields.  This functions generates strings representing genotypes in 
that order, formatted for the \code{GT} field of a VCF.
}
\usage{
genotypeStrings(ploidy, n_alleles, sep = "/")
}
\arguments{
  \item{ploidy}{
An integer indicating the ploidy.
}
  \item{n_alleles}{
An integer indicating the number of alleles, including the reference allele.
}
  \item{sep}{
A character for separating alleles within genotype strings.
}
}
\value{
Using the default \code{sep} value, a character vector containing unphased 
genotypes formatted for the \code{GT} field in a VCF file, in the order
indicated in the VCF specification for the \code{GL} field.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{enumerateGenotypes}}
}
\examples{
# diploid with one reference and one alternative allele
genotypeStrings(2, 2)

# tetraploid with one reference and two alternative alleles 
genotypeStrings(4, 3)
}

\keyword{ utilities }
