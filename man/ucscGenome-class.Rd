\name{ucscGenome-class}
\Rdversion{1.1}
\docType{class}
\alias{ucscGenome-class}
\alias{ucscGenome}
\alias{addEnsembl}
\alias{addEnsembl-methods}
\alias{addEnsembl,ucscGenome-method}
\alias{addIsoforms}
\alias{addIsoforms-methods}
\alias{addIsoforms,ucscGenome-method}
\alias{addXref}
\alias{addXref-methods}
\alias{addXref,ucscGenome-method}
\alias{extractTranscript,ucscGenome-method}
\alias{getXref}
\alias{getXref-methods}
\alias{getXref,ucscGenome-method}
\alias{tableFeatures,ucscGenome-method}
\alias{tableTranscript.id,ucscGenome-method}

\title{Class \code{"ucscGenome"}}
\description{
ucscGenome class: Represents data stored for UCSC genome. The standard way to import data is to download a "gtf" file from the UCSC Genome Browser (-> Table Browser). Download the "knownGene" Table in output format "GTF". Then import the data via the \code{read.gtf} function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{ucscGenome()}. 
}
\section{Slots}{
  \describe{
    \item{\code{basedir}:}{Object of class \code{"character"} Directory where SQLite database is written.  }
    \item{\code{ev}:}{Object of class \code{"environment"} Environment that contains data structures. Optionally, there are gtf, attr and additionally xref data.frames. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "refGenome")}: Creates a sensible printout. }
    \item{getGtf}{\code{signature(object = "refGenome")}: Returns content of gtf table. }
    \item{setGtf}{\code{signature(object = "refGenome")}: Writes content of gtf table. }
    \item{getAttr}{\code{signature(object = "refGenome")}: Returns content of attribute table. }
    \item{setAttr}{\code{signature(object = "refGenome")}: Writes content of attribute table. }
    \item{read.gtf}{\code{signature(object = "refGenome")}: Imports content of gtf file. This is the basic mechanism for data import. It works the same way for ucscGenome and for ensemblGenome. The attribute items are parsed and written in parts to gtf table and attribute table.  }
    \item{writeDB}{\code{signature(object = "refGenome")}: Copies content of gtf, attr and xref table to database. }
    \item{addEnsembl}{\code{signature(object = "ucscGenome")}: Imports UCSC 'knownToEnsembl' table. It's appended to the gtf table.}
    \item{addIsoforms}{\code{signature(object = "ucscGenome")}: Imports UCSC ' knownIsoforms' table. It's appendet to the gtf table.}
    \item{addXref}{\code{signature(object = "ucscGenome")}: Imports UCSC 'kgXref' table. A 'geneSymbol' column is added to gtf table. The rest is written into xref table. }
        \item{extractByGeneName}{\code{signature(object="ucscGenome",geneNames="character")}: Extracts ucscGenome object which contains table subsets. When none of the \code{geneNames} matches, the function returns \code{NULL}.}
    \item{getXref}{\code{signature(object = "ucscGenome")}: Returns content of xref table. }
        \item{getGenePositions}{\code{signature(object="ucscGenome",force="logical")}: Extracts table with position data for whole genes (smallest exon start position and largest exon end position. A copy of the table will be placed inside the internal environment. Upon subsequent call only a copy of the contained table is returned unless \code{force=TRUE} is given. Upon \code{force=TRUE} new gene positions are calculated regardless of existing tables.)}
    \item{loadGenome}{\code{signature(filename = "character")}: Imports data from stored R-Environment Image. }
    \item{loadGenomeDb}{\code{signature(filename = "character")}: Imports content of object from sqlite3 database. }
    \item{tableFeatures}{\code{signature(object="ucscGenome")}: Tables content of "feature" column.}
    \item{tableTranscript.id}{\code{signature(object="ucscGenome")}: Tables values in transcript_id column.}
    \item{extractTranscript}{\code{signature(object="ucscGenome",transcripts="character")}: Extracts an object which contains data for subset defined by transcript names.}
	 }
}
\references{http://genome.ucsc.edu/}
\author{Wolfgang Kaisers}
\examples{
# + + + + + + + + + + + + + + + + + + + #
# Loading and saving
# From and to R-image (fast loading)
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
uc
saveGenome(uc,"uc.RData",useBasedir=FALSE)
ucr<-loadGenome("uc.RData")
# + + + + + + + + + + + + + + + + + + + # 
# Extract data for Primary Assembly seqids
ucpa<-extractSeqids(uc,ucPrimAssembly())
# Extract data for indival Genes
ddx<-extractByGeneName(uc,"DDX11L1")
ddx
# Extract range limits of entire Genes
gp<-getGenePositions(uc)
gp
tableFeatures(uc)
extractByGeneName(ucpa,"DDX11L1")
tableTranscript.id(ucpa)
# + + + + + + + + + + + + + + + + + + + #
# Create object from scratch
# uc<-ucscGenome()
# basedir(uc)<-"/my/genome/basedir"
# Place all UCSC-files in folder
# read.gtf(uc,"knownGene.gtf")
# addXref(uc,"kgXref.csv")
# addEnsembl(uc,"knownToEnsembl.csv")
# addIsoforms(uc,"knownisoforms.csv")
# + + + + + + + + + + + + + + + + + + + #
}
\keyword{classes}
