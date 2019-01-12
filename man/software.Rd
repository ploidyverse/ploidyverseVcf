\name{software}
\alias{software}
\alias{software<-}
\alias{software,VCF-method}
\alias{software<-,VCF-method}
\title{
Identify Software Used for Calling Genotypes in a VCF
}
\description{
The \code{software} function retrieves information about the software used to
generate genotype calls present in a ploidyverse \code{\link[=VCF-class]{VCF}}.
The corresponding assignment function can be used to add the same data to a
\code{VCF} object.
}
\usage{
software(object)

software(object) <- value
}
\arguments{
  \item{object}{
A \code{VCF} object.
}
\item{value}{
A named character vector, named list, or data frame.  The names (or column
names) \code{Software}, \code{Version}, \code{Model}, and 
\code{Description} must be present.  The name \code{ID} must also be present
if multiple values are provided for each field.  Any additional named items
will also be retained in the output, allowing for custom fields.
}
}
\value{
\code{software} returns a \code{\link{DataFrame}} containing information about
the software used for genotype calling.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{markValidity}}
}
\examples{
\dontrun{
software(myvcf) <- c(Software = "polyRAD", Version = "1.1", Model = "IteratePopStruct",
                     Description = "Just some genotypes")
software(myvcf)
}
}

\keyword{ methods }
