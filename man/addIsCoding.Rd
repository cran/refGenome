\name{addIsCoding}
\alias{addIsCoding}
\alias{addIsCoding-methods}
\alias{addIsCoding,ensemblJunctions-method}
\title{Add information on coding/non-coding status.}
\description{The \code{addIsCoding} method extracts positions from 'CDS'
features in the \code{ensemblGenome} object.
The positions of the 'CDS' features are compared with the splice-junction
positions in the \code{ensemblJunctions} object. When a match is found,
the splice junction is marked as coding, otherwise the junction is marked as
non-coding. 
The coding information is calculated for each flanking exon of the 
splice-junction.
The column name for the left side (lend) is \code{licd}, the name for the
right side (rstart) is \code{ricd}.}
\usage{
addIsCoding(object,ens)
}
\arguments{
    \item{object}{\code{ensemblJunctions}.
        The object to which \code{isCoding} columns are added.}
    %
    \item{ens}{\code{ensemblGenome}.
        Genome object which should contain 'CDS' features.}
}
\author{Wolfgang Kaisers}
\examples{
ef<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ef)
enj<-getSpliceTable(ens)
addIsCoding(enj,ens)
}
\keyword{hbond}
