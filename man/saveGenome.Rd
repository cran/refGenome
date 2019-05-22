\name{saveGenome}
\alias{saveGenome}
\alias{saveGenome-methods}
\alias{saveGenome,refGenome-method}
\alias{loadGenome}
\title{Saving and loading \code{refGenome} objects}
\description{
\code{refGenome} objects contain all annotation data inside an environment.
During saving and loading, the whole content of the environment is loaded and
saved.
The correct object type is also determined in this way.
}
\usage{
saveGenome(object,filename,useBasedir=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{refGenome (or derived)}
    %
    \item{filename}{Character. Filename of the R-data-file wich is written.}
    %
    \item{useBasedir}{Logical. When TRUE the output file is written into
        basedir.
        Otherwise the directory depicted by \code{dirname(filename)} is used.}
    %
    \item{...}{(unused)}
}
\author{Wolfgang Kaisers}
\examples{
##-------------------------------------##
## A) Ensembl genome
##-------------------------------------##
ensfile <- system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens <- loadGenome(ensfile)
\donttest{
saveGenome(ens, "ens.RData", useBasedir=FALSE)
}

## B) Ensembl junctions
junc <- getSpliceTable(ens)
\donttest{
saveGenome(junc, "junc.RData", useBasedir=FALSE)
loadGenome("junc.RData")
}

##-------------------------------------##
## C) Ensembl exons
##-------------------------------------##
enex <- refExons(ens)
\donttest{
saveGenome(enex, "enex.RData", useBasedir=FALSE)
er <- loadGenome("enex.RData")
}

##-------------------------------------##
## D) UCSC genome:
##-------------------------------------##
ucfile <- system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc <- loadGenome(ucfile)
\donttest{
saveGenome(uc, "uc.RData", useBasedir=FALSE)
}
##-------------------------------------##
## E) UCSC junctions
##-------------------------------------##
junc <-getSpliceTable(uc)
\donttest{
saveGenome(junc, "junc.RData", useBasedir=FALSE)
jr <- loadGenome("junc.RData")
}

##-------------------------------------##
## F) UCSC exons
##-------------------------------------##
ucex <- refExons(uc)
\donttest{
saveGenome(ucex, "ucex.RData", useBasedir=FALSE)
ur <- loadGenome("ucex.RData")
}
}
\keyword{overlap}
