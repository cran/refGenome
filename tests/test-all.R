

##============================================================================##
## A    Load prerequisites
##============================================================================##

require(refGenome)


##============================================================================##
## B    Initialize example data
##============================================================================##

ensfile <- system.file("extdata",
                       "hs.ensembl.62.small.RData",
                       package = "refGenome", mustWork=TRUE)

##----------------------------------------------------------------------------##
## B.1 Load Ensembl genome
##----------------------------------------------------------------------------##
ens <- loadGenome(ensfile)
enex<-refExons(ens)
gp <- getGenePositions(ens)
junc<-getSpliceTable(ens)



##============================================================================##
## C    Run tests
##============================================================================##

##----------------------------------------------------------------------------##
## C.1 Test overlap juncs
##----------------------------------------------------------------------------##

##                                                                            ##
## Requires: Initialized objects (as done by test-all.R header)
##                                                                            ##

##                                                                            ##
## C.1.1 Overlap juncs
##                                                                            ##

##                      1       2       3       4       5       6       7 ##
qry<-data.frame(id = 1:7, seqid = "1",
                lstart = c(10100L, 11800L, 12220L, 12220L, 12220L, 32000L, 40000L),
                lend =   c(10100L, 12000L, 12225L, 12227L, 12227L, 32100L, 40100L),
                rstart = c(10200L, 12200L, 12057L, 12613L, 12650L, 32200L, 40200L),
                rend =   c(10300L, 12250L, 12179L, 12620L, 12700L, 32300L, 40300L))
##                      1       2       3       4       5       6       7 ##

##                                                                            ##
## C.1.2
##                                                                            ##
res<-overlapJuncs(qry,junc)

if(! all(is.na(res$sod)==c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)) )
    stop("[test_overlap_juncs] Wrong res$sod NA's.")

if(sum(res$nref) != 27)
    stop("[test_overlap_juncs] Wrong sum of res$nref.")


##                                                                            ##
## C.1.3    Test overhang codes (No requirements)
##                                                                            ##
rid     <- 1L
rLstart <- 1000L
rLend   <- 2000L
rRstart <- 3000L
rRend   <- 4000L
rMaxEnd <- 4000L
#
qLstart <- as.integer(c(0700, 0800, 0900, 1000, 1000, 1000, 1000))
qLend   <- as.integer(c(0800, 0900, 2000, 2000, 2000, 2000, 5000))
qRstart <- as.integer(c(0900, 3000, 3000, 3000, 3000, 5000, 6000))
qRend   <- as.integer(c(4000, 4000, 4000, 4000, 5000, 7000, 7000))
qid     <- as.integer(1:length(qLstart))
dfr <- .Call(refGenome:::C_gap_overlap, qid, qLstart, qLend, qRstart, qRend,
             rid, rLstart, rLend, rRstart, rRend, rMaxEnd)

if(!all(dfr$ovhl==c("inp","int","ext",rep("no",4))))
    stop("[test_overlap_juncs] Wrong ovhl codes")
if(!all(dfr$ovhr==c(rep("no",4),"ext","int","inp")))
    stop("[test_overlap_juncs] Wrong ovhr codes")



##----------------------------------------------------------------------------##
## C.2 Test read.gtf features
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## hs.ensembl.76.small.gtf
## contains the first 100 lines of file
## ftp://ftp.ensembl.org/pub/release-76/gtf/
##                                  homo_sapiens/Homo_sapiens.GRCh38.76.gtf.gz
##----------------------------------------------------------------------------##

ens76 <- system.file("extdata",
                     "hs.ensembl.76.small.gtf", package = "refGenome", mustWork=TRUE)

en76 <- ensemblGenome()
basedir(en76) <- dirname(ens76)
read.gtf(en76, basename(ens76))

if(!exists("genes", where=en76@ev, inherits=FALSE))
    stop("[test_read_gtf] 'genes' table does not exist")

if(!all(getGeneTable(en76)$gene_name[1:2]==c("DDX11L1","WASH7P")))
    stop("[test_read_gtf] Wrong gene names in 'genes' table.")


##----------------------------------------------------------------------------##
## C.3 Test unify ranges
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Output consistency of unify Ranges:
## Consecutive ranges must not ovelap.
##----------------------------------------------------------------------------##

urg <- unifyRanges(enex)
clevels <- levels(urg@ev$gtf$seqid)
n <- length(clevels)

for(i in 1:n)
{
    xtr <- urg@ev$gtf[urg@ev$gtf$seqid == clevels[i], ]
    xb <- c(xtr$begin, 0)
    xe <- c(0,xtr$end)
    diff <- xb - xe
    if(table(diff > 0)[1] > 1)
        stop("[unifyRanges] inconsistent output!")
}

##----------------------------------------------------------------------------##
## C.4 Test unify juncs
##----------------------------------------------------------------------------##


##----------------------------------------------------------------------------##
# Artificially added transcript_biotype data from ensembl version 80
# For transcript_id's which are not found idn ensembl 80,
# the transcript_biotype is manually set to "nonsense_mediated_decay".
##----------------------------------------------------------------------------##
# ensfile <- system.file("extdata", "hs.ensembl.62.small.RData",
#                                         package = "refGenome", mustWork=TRUE)
# ens <- loadGenome(ensfile)
# jtn <- ens@ev$gtf$transcript_id
# mtc <- match(jtn, en80@ev$gtf$transcript_id)
# biot <- en80@ev$gtf$transcript_biotype[mtc]
# biot[is.na(biot)] <- "nonsense_mediated_decay"
# ens@ev$gtf$transcript_biotype <- biot
# saveGenome(ens,
#     file="~/projects/R/refGenome/inst/extdata/hs.ensembl.62.small.RData",
#     useBasedir=FALSE)
#
##----------------------------------------------------------------------------##


uj <- unifyJuncs(junc)

if(! all(uj@ev$gtf$cnNmd >= 0) )
    stop("[test_unify_juncs] All cnNmd must be >= 0")

if(! all(uj@ev$gtf$cnNmd <= uj@ev$gtf$nSites) )
    stop("[test_unify_juncs] All cnNmd must be <= nSites")

rm(uj)


##============================================================================##
## D    Cleanup
##============================================================================##


rm(ensfile, ens, gp, junc)
gc()

cat("[refGenome] rest-all.R tests finished.\n")


##============================================================================##
## END OF FILE
##============================================================================##
