\name{nGen}
\alias{alleleCopy}
\alias{dmultinom}
\alias{dDirichletMultinom}
\alias{enumerateGenotypes}
\alias{genotypeFromIndex}
\alias{indexGenotype}
\alias{makeGametes}
\alias{nGen}
\alias{selfingMatrix}
\title{
Utilities for Multiallelic Genotypes in Rcpp
}
\description{
These functions provide some utilities for processing multiallelic genotypes,
including multinomial probability, adjustment of genotype frequencies under
certain mating schemes, and ordering of genotypes according to the 
\href{https://github.com/samtools/hts-specs}{VCF specification}.  They are 
all implemented in Rcpp, so that they can be used from R or from Rcpp.
}
\usage{
dmultinom(x, prob)
dDirichletMultinom(x, prob, alpha)

nGen(ploidy, nalleles)
enumerateGenotypes(ploidy, nalleles)
indexGenotype(genotype)
genotypeFromIndex(index, ploidy)
alleleCopy(genotype, nalleles)

makeGametes(genotype)
selfingMatrix(ploidy, nalleles)
}
\arguments{
\item{x}{
A vector of integers indicating number of items of each category sampled,
for example copy numbers of alleles within a genotype (when estimating
genotype frequencies under HWE), or sequence read depth for each allele
(when estimating genotype likelihood).
}
\item{prob}{
A numeric vector the same length as \code{x}, indicating the probability
of sampling an item from each category.
}
\item{alpha}{
An overdispersion parameter.  Lower values result in more overdispersed
results, whereas higher values will cause the output of 
\code{dDirchletMultinom} to more closely resemble that of \code{dmultinom}.
}
  \item{ploidy}{
An integer indicating the ploidy of the organism.
}
  \item{nalleles}{
An integer indicating the number of alleles, including the reference allele.
}
\item{genotype}{
An integer vector indicating a genotype.  The length of the vector should be
the same as the ploidy of the organism.  Each item in the vector represents
one copy of the locus, with the integer representing an allele.  The reference
allele should be represented as zero, and alternative alleles as integers above
zero.  The vector must be sorted in ascending order.
}
\item{index}{
The index (starting at zero) of a genotype in the VCF order.
}
}
\details{
Since these functions are designed to be run repeatedly in computationally
intensive loops, they don't perform error checking.  It is highly advisable
to perform error checking at a higher level in the code, for example checking
once that the ploidy is above zero before making millions of calls to 
\code{nGen} with that ploidy value.

\code{dmultinom} should give identical results to 
\code{\link[stats:Multinom]{stats::dmultinom}}, but be faster and accessible
from within Rcpp.

Where \eqn{K} is the length of \eqn{x}, \eqn{p} is the vector of 
probabilities, \eqn{\alpha} is the overdispersion parameter, and
\eqn{n} is the sum of \eqn{x}:

\deqn{n = \sum_{i = 1}^K x_i}{n = sum(i = 1, ..., K) x[i]}

The output of \code{dmultinom} is

\deqn{
n! \prod_{i = 1}^K \frac{p_i^{x_i}}{x_i!}}{
n! * prod(i = 1, ..., K)((p[i] ^ x[i])/x[i]!)
}

And the output of \code{dDirichletMultinom} is

\deqn{
\frac{(n!) \Gamma(\alpha)}{\Gamma(n + \alpha)} 
\prod_{i = 1}^K \frac{\Gamma(p_i * \alpha + x_i)}{(x_i!) \Gamma(p_i * \alpha)}}{
n! * \Gamma(\alpha) / \Gamma(n + \alpha) * 
prod(i = 1, ...K)(\Gamma(p[i] * \alpha + x[i])/ x[i]! / \Gamma(p[i] * \alpha))
}
}
\value{
\code{dimultinom} and \code{dDirchletMultinom} return a single numeric value
between 0 and 1 indicating the multinomial or Dirichlet-multinomial 
probability, repectively.

\code{nGen} returns an integer indicating how many unique genotypes are 
possible for a given ploidy and number of alleles.

\code{enumerateGenotypes} returns an integer matrix, with as many columns as
the ploidy, listing all possible genotypes.  Each row represents one genotype.
Each element represents an allele in the genotype, with zero indicating the
reference allele and higher integers indicating alternative alleles.  The
genotypes are ordered according to the VCF specification.

\code{indexGenotype} returns an integer indicating the index of a genotype
within the VCF specification (i.e. in what row of the matrix returned by
\code{enumerateGenotypes} would we find that genotype), starting at index zero.

\code{genotypeFromIndex} returns an integer vector indicating the genotype
at a given index (starting at index zero, and according to the VCF
specification).

\code{alleleCopy} returns an integer vector indicating the genotype in an
alternative format: there is one element for each allele (including the
reference allele and all alternative alleles), and the value of the element
indicates how many copies of that allele that genotype has.  This format
can be useful for input to \code{dmultinom}.

\code{makeGametes} returns an integer matrix listing all possible gametes.
The number of columns will be equal to half the ploidy of the input genotype.
Each row represents one gamete.  Gametes may be listed multiple times in
order to reflect their relative frequency.  Fully polysomic inheritance is
assumed with no double reduction.

\code{selfingMatrix} returns a numeric matrix with parent genotypes in rows and
progeny genotypes in columns, in VCF order.  Matrix elements indicate the
frequency with which each parent genotype will produce each progeny genotype
after self-fertilization.  Polysomic inheritance and no double reduction are
assumed.  If \eqn{f} is a vector of genotype frequencies and \eqn{A} is
the matrix output by this function, \eqn{fA} is the vector of genotype
frequencies after one round of self-fertilization.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{genotypeStrings}}
}
\examples{
# Likelihood of genotype 0/1/2/2 having 20 reads for 0, 25 reads
# for 1, and 35 reads for 2.
dmultinom(c(20, 25, 35), c(0.25, 0.25, 0.5))
dDirichletMultinom(c(20, 25, 35), c(0.25, 0.25, 0.5), 9)

# A tetraploid with two alternative alleles
nGen(4, 3)
enumerateGenotypes(4, 3)
# A diploid with one alternative allele
nGen(2, 2)
enumerateGenotypes(2, 2)

# Indices for two diploid genotypes
indexGenotype(c(0, 2))
indexGenotype(c(1, 1))
# Indices for two tetraploid genotypes
indexGenotype(c(0, 0, 0, 0))
indexGenotype(c(0, 0, 2, 2))

# Convert a tetraploid genotype to copy number format
alleleCopy(c(1, 1, 1, 2), 4)

# Generate gametes
makeGametes(c(0, 1)) # Diploid heterozygous
makeGametes(0:3) # Tetraploid fully heterozygous
makeGametes(0:5) # Hexaploid fully heterozygous
makeGametes(c(0,0,0,1)) # Tetraploid partial heterozygote

# Self-fertilization matrix for a biallelic diploid
selfingMatrix(2, 2)
}
\keyword{ utilities }
