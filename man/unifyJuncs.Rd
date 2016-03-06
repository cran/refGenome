\name{unifyJuncs}
\Rdversion{1.1}
\alias{unifyJuncs}
\alias{unifyJuncs-methods}
\alias{unifyJuncs,refJunctions-method}
\title{Unification of splice-junctions.}
\description{
Overlaps query ranges with reference ranges. The function assumes that
there is no overlap between reference ranges.
}
\usage{
unifyJuncs(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{refJunctions} object. Contains splice-junction data.}
}
\details{
Many splices sites are multiple times contained when calculated from 
transcripts.
In order to obtain unique splice positions the function extracts one data
set per (seqid, lend, rstart) combination.
For each site the annotation information (gene_id, strand, fexid)
is extracted from the most abundand gene name.
\code{cnNmd} suffices: 0 <= cnNmd <= nSites. When cnNmd=0, the junction
is only present in Transcripts with biotype \code{nonsense mediated decay}.
}
\value{
The returned object has the same class as the passed argument (refJunctions
or derived).
Therefore, it's possible to use unified junctions in the same
way in downstream analysis as un-unified junctions.
The contained gtf - data.frame contains the following columns:
\tabular{ll}{
    id \tab Numeric index for unique site.\cr
    seqid \tab Chromosome identifier\cr
    lstart \tab left start\cr
    lend \tab left end\cr
    rstart \tab right start\cr
    rend \tab right end \cr
    nSites \tab Number of \code{refJunctions} (transcripts) 
        that contain this site.\cr
    %
    gene_id \tab Gene identifier.\cr
    strand \tab Strand on which gene resides.\cr
    fexid \tab Id value of first \code{refJunction} with coordinates.\cr
    cnNmd \tab Number of \code{refJunctions} whose transcript_biotype is not
        \code{nonsense_mediated_decay}.\cr
}
}
\author{
Wolfgang Kaisers
}
\examples{
# + + + + + + + + + + + + + + + + + + + #
# A) Ensembl
# + + + + + + + + + + + + + + + + + + + #
ef<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ef)
enj<-getSpliceTable(ens)
ufe<-unifyJuncs(enj)
saveGenome(enj,"enj.RData",useBasedir=FALSE)
enjr<-loadGenome("enj.RData")
# + + + + + + + + + + + + + + + + + + + #
# B) UCSC
# + + + + + + + + + + + + + + + + + + + #
uf<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(uf)
ucj<-getSpliceTable(uc)
ufu<-unifyJuncs(ucj)
saveGenome(ucj,"ucj.RData",useBasedir=FALSE)
ucjr<-loadGenome("ucj.RData")
}
\keyword{refJunctions}
