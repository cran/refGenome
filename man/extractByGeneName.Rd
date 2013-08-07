\name{extractByGeneName}
\alias{extractByGeneName}
\alias{extractByGeneName-methods}
\alias{extractByGeneName,refGenome-method}
\title{Extract subsets of refGenome by gene-name.}
\description{
The function takes objects derived from \code{refGenome} and returns a subset in which \code{gene_name} matches the given values. The returned object is of the same class as the given object.
}
\usage{
extractByGeneName(object,geneNames,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{refGenome (or derived). Object from which subset is extracted.}
  \item{geneNames}{Character. Gene names }
  \item{...}{(unused)}
}
\value{Same class as \code{object}}
\author{Wolfgang Kaisers}
\examples{
ensfile<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)
ws7<-extractByGeneName(ens,"WASH7P")
ws7
}
\keyword{overlap}
