\name{array3D_to_matrixList}
\alias{array3D_to_matrixList}
\alias{matrixList_to_array3D}
\title{
Convert Between Two-Dimensional Lists and Three-Dimensional Arrays
}
\description{
The \code{\link[=VCF-class]{VCF}} class from \pkg{VariantAnnotation} stores
genotype data in matrices, with one matrix for each field.  Fields such as
\code{AD}, \code{GP}, and \code{GL} can store multiple values, in which case
the matrix is a two-dimensional list of vectors.  Where all vectors are the
same length, however, it can be efficient to work with the data as a 
three-dimensional array in R.

These two functions convert genotype data between these two formats, to
make it easier for users and package developers to get their data into and
out of \code{VCF} objects.
}
\usage{
array3D_to_matrixList(arr)

matrixList_to_array3D(mat)
}
\arguments{
  \item{arr}{
A three-dimensional array with samples in the second dimension and loci in the
third dimension.  The first dimension can represent unique genotypes, alleles, 
haplotypes, or allele copy number.
}
  \item{mat}{
A matrix with loci as rows and samples as columns.  Each item in the matrix is
a vector.  All vectors must be the same length.
}
}
\details{
Sample and locus names are carried over from input to output (\code{dimnames} 
of \code{arr} and \code{mat}).
}
\value{
\code{array3D_to_matrixList} returns a matrix formatted like \code{mat}.

\code{matrixList_to_array3D} returns an array formatted like \code{arr}.  The
names of the first dimension will be integers converted to character, starting
with zero.
}

\author{
Lindsay V. Clark
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# create a dummy dataset for this example
dummy <- array(rnorm(30), dim = c(5, 2, 3),
               dimnames = list(as.character(0:4), c("sam1", "sam2"),
                               c("loc1", "loc2", "loc3")))

# convert to matrix
dummy_mat <- array3D_to_matrixList(dummy)

# example of how data are stored differently
dummy[,1,2]
dummy_mat[[2,1]]

# convert to array
dummy_arr <- matrixList_to_array3D(dummy_mat)

identical(dummy, dummy_arr) # TRUE
}

\keyword{ manip }
