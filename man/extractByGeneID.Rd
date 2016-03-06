\name{extractByGeneId}
\alias{extractByGeneId}
\alias{extractByGeneId-methods}
\alias{extractByGeneId,refGenome,character-method}
\alias{extractByGeneId,refJunctions,character-method}
\title{Extract subsets of refGenome by gene-ID.}
\description{
The function takes objects derived from \code{refGenome} or \code{refJunctions}
and returns a subset in which \code{gene_name} matches the given values.
The returned object is of the same class as the given object.
}
\usage{
extractByGeneId(object,geneids,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{refGenome (or derived).
        Object from which subset is extracted.}
    %
    \item{geneids}{Character. Vector with gene ID's.}
    %
    \item{...}{(unused)}
}
\value{Same class as \code{object}}
\author{Wolfgang Kaisers}
\examples{
# + + + + + + + + + + + + + + + + + + #
# A) Extract from Genome
# + + + + + + + + + + + + + + + + + + #
ensfile<-system.file("extdata", 
                    "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)
ddx<-extractByGeneId(ens,"ENSG00000223972")
ddx
# + + + + + + + + + + + + + + + + + + #
# B) Extract from splice junctions
# + + + + + + + + + + + + + + + + + + #
junc<-getSpliceTable(ens)
ddx<-extractByGeneId(junc,"ENSG00000223972")
ddx
}
\keyword{overlap}
