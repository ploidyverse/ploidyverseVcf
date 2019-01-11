\name{sampleinfo}
\alias{sampleinfo}
\alias{sampleinfo<-}
\alias{sampleinfo,VCF-method}
\alias{sampleinfo<-,VCF-method}
\title{
Access or Modify Sample Metadata in a VCF
}
\description{
The \code{sampleinfo} function retrieves a \code{\link{DataFrame}} containing
sample metadata from a \code{\link[=VCF-class]{VCF}} object, stored in
\code{meta(header(object))$SAMPLE}, or the \code{##SAMPLE} lines of the VCF
file.  The corresponding assignment function can be used to add this same
metadata to a \code{VCF} object.
}
\usage{
sampleinfo(object)

sampleinfo(object) <- value
}

\arguments{
  \item{object}{
A \code{VCF} object.
}
\item{value}{
A \code{data.frame}, \code{DataFrame}, or other object that can be coerced to a
\code{DataFrame}.  Row names should correspond to 
\code{samples(header(object))}, the names of the samples.  Columns contain
metadata about each sample.
}
}
\details{
\code{Species} and \code{Ploidy} are expected to be present in the column names
of \code{value}; a warning will be given if they are not.  The ploidy column 
should be formatted \dQuote{2x} for a diploid, \dQuote{4x} for an 
autotetraploid, \dQuote{2x+2x} for an allotetraploid, etc.  Ploidy should be 
with respect to the reference genome if used, or with respect to the monoploid
genome if data were generated without a reference sequence.

Other columns to optionally add are \code{CollectionLocation}, 
\code{MaterialProvider}, \code{PopulationDesign}, \code{BioSample}, and 
\code{DOI}.  These metadata fields are described in the \dQuote{Variant Call
Format specification for the ploidyverse} document in the \dQuote{doc} folder
of this package installation.  For any of these, as well as \code{Species}
and \code{Ploidy}, descriptions are automatically added to 
\code{meta(header(object))$META} (corresponding to \code{##META} lines in the
file).

Custom columns can be included, in which case a warning will be printed
instructing the user to add column descriptions.  If your package uses custom
column names and you want to avoid the warning, let Lindsay Clark know, and
add a commit or pull request to add your columns to the list in the method
definition.

If \code{value} is a \code{data.frame}, it should typically be constructed with
the argument \code{stringsAsFactors = FALSE}.
}
\value{
\code{sampleinfo} returns a \code{DataFrame} containing sample metadata.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{markValidity}}
}
\examples{
# To be added
}

\keyword{ methods }

