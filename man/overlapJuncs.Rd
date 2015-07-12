\name{overlapJuncs}
\alias{overlapJuncs}
\title{overlapJuncs function}
\description{
Overlaps query gap-sites (from BAM alignment data) with annotated 
splice junctions (from reference genome annotation).
}
\usage{
overlapJuncs(qry, junc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{qry}{\code{data.frame}.
        Table with query ranges. 
        qry should have columns 'id', 'seqid', 'lstart', 
        'lend', 'rstart', 'rend'.}
    %
    \item{junc}{\code{refJunctions}. 
        Object which contains table of splice junctions in reference genome.}
    %
}
\details{
The function finds optimal overlapping hits for alignment gap-sites in
annotated splice-sites. 
A gap-site is the combination of two genomic regions (= exons) which enclose
an intermediary (= intron).
The function identifies junction records which overlap with the given
gap-site (=hits) and select a junction with an optimal fit.
The goodness of fit is measured by the distance of the inner
gap boundaries (= the splice sites) between query and junction record.
A junction with minimal sum of upstream and downstrem distances is selected.
The selection of the best hit depends on the order a version of the junction
table wich is sorted by \code{lstart} and \code{rend}.)
}
\value{
The function returns a data.frame
    \item{qid}{\code{Integer}. Query id value from qry table}
    %
    \item{refid}{\code{Integer}. Reference id from junctions object for best hit.}
    %
    \item{ldiff}{\code{Integer}. Difference between lend values in qry and 
        \code{junc} table for best hit (refid) record.}
    %
    \item{rdiff}{\code{Integer}. Difference between rstart values in qry
        and \code{junc} table for best hit (refid) record.}
    %
    \item{nref}{\code{Integer}. Number of \code{junc} records which possibly
        overlap with query item. \code{nref=0} when no overlap
        has been found for query.}
    %
    \item{sod}{\code{Integer}. Sum of distances (=abs(ldiff) + abs(rdiff)).
        \code{sod=0} when qry exactly hits an annotated site. \code{sod=NA}
        when no overlap has been found for query.}
    %
    \item{first_refid}{\code{Integer}. 
        \code{id} for first overlapping record in \code{junc} table.}
    %
    \item{last_refid}{\code{Integer}. 
        \code{id} for last overlapping record in \code{junc} table.}
    %
    \item{nadv}{\code{Integer}. 
        Number of advancing iterations during search for }
    %
    \item{strand}{Strand value derived from annotation.}
    %
    \item{gene_id}{Gene id from \code{refJunctions} or \code{genpos} table.}
    %
    \item{transcript_id}{Transcript id from \code{refJunctions} table.}
    %
    \item{gene_name}{Gene name from \code{refJunctions} or \code{genpos} table.}
    %
}
\author{
Wolfgang Kaisers
}
\examples{
# + + + + + + + + + + + + + + + + + + #
# A) Example query  data
# + + + + + + + + + + + + + + + + + + #
##                          1       2       3       4       5       6       7 ##
qry <- data.frame(id = 1:7, seqid = "1",
            lstart = c(10100L, 11800L, 12220L, 12220L, 12220L, 32000L, 40000L),
            lend =   c(10100L, 12000L, 12225L, 12227L, 12227L, 32100L, 40100L),
            rstart = c(10200L, 12200L, 12057L, 12613L, 12650L, 32200L, 40200L),
            rend =   c(10300L, 12250L, 12179L, 12620L, 12700L, 32300L, 40300L))
##                          1       2       3       4       5       6       7 ##

# + + + + + + + + + + + + + + + + + + #
# B) Example reference data
# + + + + + + + + + + + + + + + + + + #
# B.1) Ensembl genome:
ensfile <- system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens <- loadGenome(ensfile)
gp <- getGenePositions(ens)
# B.2) Ensembl junctions:
junc <- getSpliceTable(ens)
# + + + + + + + + + + + + + + + + + + #
# C) Do overlap
# + + + + + + + + + + + + + + + + + + #
res <- overlapJuncs(qry, junc)
}

\keyword{overlap}
