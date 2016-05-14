\name{getGeneId}
%
\alias{getGeneId}
\alias{getGeneId-methods}
\alias{getGeneId,ensemblGenome,character-method}
\alias{getGeneId,ensemblGenome,factor-method}
%
\title{Extract gene ID from gene-name and getGeneId.}
\description{
The function searches in the given refGenome object for \code{gene_name} values
and returns gene id values. For unkown gene names, NA is returned.
}
\usage{
getGeneId(object, geneName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{ensemblGenome,
        Object from which gene ID are extracted.}
    %
    \item{geneName}{Character: Vector of gene names (factor also accepted).}
}
\value{character}
\author{Wolfgang Kaisers}
\examples{
# + + + + + + + + + + + + + + + + + + #
# A) Ensembl genome:
# + + + + + + + + + + + + + + + + + + #
ensfile<-system.file("extdata", 
                        "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)
getGeneId(ens, c("DDX11L1", "WASH7P"))
}
\keyword{overlap}
