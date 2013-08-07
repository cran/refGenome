\name{writeDB}
\alias{writeDB}
\alias{writeDB-methods}
\alias{writeDB,refGenome-method}
\alias{loadGenomeDb}
\title{Saving and loading \code{refGenome} objects to and from SQLite databases.}
\description{
\code{refGenome} objects keep annotation data in data.frames. The content of the data.frames is written to or extracted from a SQLite database.
}
\usage{
writeDB(object,filename,useBasedir=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{refGenome (or derived)}
  \item{filename}{Character. Filename of the R-data-file wich is written.}
  \item{useBasedir}{Logical. When TRUE the database will be created in basedir. Otherwise the database will be created in the current working directory or elsewhere (when path is included in filename).}
  \item{...}{(unused)}
}
\author{Wolfgang Kaisers}
\examples{
# + + + + + + + + + + + + + + + + + + + #
# A) Ensembl genome:
ensfile<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)
writeDB(ens,"ens.db3",useBasedir=FALSE)
ens<-loadGenomeDb("ens.db3")
# + + + + + + + + + + + + + + + + + + + #
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
writeDB(uc,"uc.db3",useBasedir=FALSE)
uc<-loadGenomeDb("uc.db3")
# + + + + + + + + + + + + + + + + + + + #
}
\keyword{overlap}
