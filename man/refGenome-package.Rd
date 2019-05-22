\name{refGenome-package}
\alias{refGenome-package}
\docType{package}
\title{
Managing annotation data for reference Genomes from UCSC and Ensembl.
}
\description{
The package contains classes for managing (GTF-) annotation data for UCSC
and Ensembl genomes.
Data can be imported, merged, viewed, searched and saved (as .RData and
as SQLite database).
There is also a C-routine for detection of overlapping (alignment) ranges
with annotated regions.
}
\details{
\tabular{ll}{
Package: \tab refGenome\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-10-06\cr
License: \tab What license is it under?\cr
Depends: \tab methods\cr
}
}
\author{
Wolfgang Kaisers
Maintainer: Wolfgang Kaisers <kaisers@med.uni-duesseldorf.de>
}
\keyword{package}
\examples{
ens<-ensemblGenome()
basedir(ens) <- system.file("extdata",package="refGenome")
ens_gtf <- "hs.ensembl.62.small.gtf"
read.gtf(ens,ens_gtf)
ddx <- extractByGeneName(ens,"DDX11L1")
ddx
fam <- extractTranscript(ens,"ENST00000417324")
fam
enpa <- extractSeqids(ens,ensPrimAssembly())
enpa
tableTranscript.id(ens)
tableTranscript.name(ens)
}
