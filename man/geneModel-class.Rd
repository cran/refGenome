\name{geneModel-class}
\Rdversion{1.1}
\docType{class}
\alias{geneModel-class}
\alias{initialize,geneModel-method}
\alias{+,geneModel,geneModel-method}
%
\alias{geneModel}
\alias{geneModel-methods}
\alias{geneModel,ensemblGenome,character-method}
\alias{geneModel,ucscGenome,character-method}
%
\alias{geneId<-}
\alias{geneId<--methods}
\alias{geneId<-,geneModel,character-method}
\alias{geneId<-,geneModel,factor-method}
%
\alias{geneName<-}
\alias{geneName<--methods}
\alias{geneName<-,geneModel,character-method}
\alias{geneName<-,geneModel,factor-method}
%
\alias{plot.geneModel}
%
\title{Class \code{"geneModel"}}
\description{
geneModel represents data for a single gene.

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{.geneModel()}.
}

\section{Slots}{
    \describe{
        \item{\code{gene_id}:}{Object of class \code{"character"}.
            ID of one single gene.}
        %
        \item{\code{gene_name}:}{Object of class \code{"character"}.
            Name of one single gene.}
        %
        \item{\code{seq_name}:}{Object of class \code{"character"}.
            Name of sequence
          (chromosome).}
        %
        \item{\code{strand}:}{Object of class \code{"character"}.
            Strand of one gene.}
        %
        \item{\code{transcripts}:}{Object of class \code{"character"}.
            Vector of transcript names.}
        %
        \item{\code{coords}:}{Object of class \code{"character"}.
            Limiting coordinates (start, end) of gene.}
        %
        \item{\code{ev}:}{Object of class \code{"environment"}
            Environment that contains data structures.
            Eventually contains list with \code{transcriptModel}
            objects.}
        %
    }
}
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "geneModel")}:
            Creates a sensible printout. }
        %
        \item{geneId}{\code{signature(object= "geneModel"}
                Returns (or sets) geneId (e.g. Ensembl ID).}
        %
        \item{geneName}{\code{signature(object= "geneModel"}
                Returns (or sets) gene name (e.g. HGNC name).}
        %
        \item{getTranscript}{\code{signature(object= "geneModel",i="numeric"}
                Returns i-th transcript object.}
        %
        \item{getTranscript}{\code{signature(object= "geneModel",i="character"}
            Returns transcript with given name.}
        %
    }
}
\usage{
geneModel(object, gene_id, interior=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{ensemblGenome or ucscGenome.
        Genome object from which gene data is extracted.}
    %
    \item{gene_id}{One single gene_id}
    %
    \item{interior}{When \code{FALSE}, exon and transcript data
        skipped.}
}
\author{
Wolfgang Kaisers
}
\examples{
# + + + + + + + + + + + + + + + + + + #
# A) Ensembl
# + + + + + + + + + + + + + + + + + + #
ensfile<-system.file("extdata",
                        "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)

gt <- getGeneTable(ens)
gene_id <- as.character(gt$gene_id[1])
gm <- geneModel(ens, gene_id, interior=TRUE)

# + + + + + + + + + + + + + + + + + + #
# B) UCSC
# + + + + + + + + + + + + + + + + + + #
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
gt <- getGeneTable(uc)
gene_id <- as.character(gt$gene_id[1])
gm <- geneModel(uc, gene_id)
# Use "+" operator
res <- gm + gm
}
\keyword{classes}
