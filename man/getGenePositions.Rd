\name{getGenePositions}
\alias{getGenePositions}
\alias{getGenePositions-methods}
\alias{getGenePositions,ensemblGenome-method}
\alias{getGenePositions,ucscGenome-method}
\alias{getGenePositions,refJunctions-method}
\title{Extract subsets of refGenome by gene-name.}
\description{
The function searches in the given data for unique \code{gene_id} values. For each item, values like gene_id, seqid and strand are extracted. Optionally (when present) also gene_name and gene_biotype are extracted. The function assignes unique id values which are ordered by gene_id values (not genetic positions). Owing to this, id values equal \code{as.numeric(gene_id)}.
}
\usage{
getGenePositions(object,by,force=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ensemblGenome, ucscGenome or ensemblJunctions. Object from which gene positions are extracted.}
  \item{by}{Character. Determines criterion by which genes are discerned. Accepted values: "gene_id" and "gene_name". For Ensembl genomes the default  is 'gene_id' and for USCC genomes the default is 'gene_name'.}
  \item{force}{Logical. When FALSE, gene positions will only be calculated when a position table is not present in local environment. The function then returns a copy.}
  \item{...}{Unused.}
}
\value{data.frame}
\details{The function stores a copy of the result in the internal environment (\code{genes} table). Upon subsequent calls the values only are re-calculated when force=TRUE is given. Otherwise the function returns a copy of the contained table. Present \code{genes} tables will be automatically saved and restored by \code{saveGenome} and \code{load.X} functions.}
\author{Wolfgang Kaisers}
\examples{
# A) Ensembl genome:
ensfile<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)
gp<-getGenePositions(ens)
# B) Ensembl junctions:
junc<-getSpliceTable(ens)
genes<-getGenePositions(junc)
# C) Ucsc genome:
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
junc<-getSpliceTable(uc)
gp<-getGenePositions(junc)
}
\keyword{overlap}
