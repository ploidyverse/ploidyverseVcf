\name{enumerateGenotypes}
\alias{enumerateGenotypes}
\alias{genotypeStrings}
\title{
Generate Genotypes Ordered According to VCF Specification
}
\description{
For any given ploidy and number of alleles, the 
\href{https://github.com/samtools/hts-specs}{VCF specification} indicates the
order in which genotypes should be listed in the \code{GP}, \code{GL}, and
\code{PL} fields.  These two functions generate genotypes in that order, as a
matrix of allele dosage and a vector of character strings, respectively.
Because these functions enumerate the set of all possible multiallelic 
genotypes, they may be useful for likelihood and probability estimates, as well
as preparing data for VCF format.
}
\usage{
enumerateGenotypes(ploidy, n_alleles)

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
\code{enumerateGenotypes} returns a matrix, with alleles represented in columns
and unique genotypes represented in rows.  Each element indicates the copy
number of that allele for that genotype.

\code{genotypeStrings} returns, by default, unphased genotypes formatted for
the \code{GT} field in a VCF file, in the same order as the rows of the matrix
returned by \code{enumerateGenotypes}.
}

\author{
Lindsay V. Clark
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# diploid with one reference and one alternative allele
enumerateGenotypes(2, 2)
genotypeStrings(2, 2)

# tetraploid with one reference and two alternative alleles 
enumerateGenotypes(4, 3)
genotypeStrings(4, 3)
}
