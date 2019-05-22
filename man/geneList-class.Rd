\name{geneList-class}
\Rdversion{1.1}
\docType{class}
\alias{geneList-class}
\alias{initialize,geneList-method}
%
\alias{geneList}
\alias{geneList-methods}
\alias{geneList,ensemblGenome,character-method}
\alias{geneList,ensemblGenome,factor-method}
\alias{geneList,ucscGenome,character-method}
\alias{geneList,ucscGenome,factor-method}
%
\alias{length,geneList-method}
\alias{show,geneList-method}
\alias{names,geneList-method}
\alias{+,geneList,geneList-method}
%
\alias{[,geneList,ANY,ANY,ANY-method}
\alias{names<-,geneList,character-method}
\alias{names<-,geneList,numeric-method}
%
\title{Class \code{"geneList"}}
\description{
geneList represents data for a single gene.

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{.geneList()}.
}
%
\section{Slots}{
    \describe{
        \item{\code{l}:}{Object of class \code{"list"}.
            List containing \code{geneModel} objects.}
    }
}
%
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "geneList")}:
            Creates a sensible printout. }
        %
        \item{length}{\code{signature(object= "geneList"}
                Returns length of contained (geneModel) list.}
        %
        \item{names}{\code{signature(object= "geneList"}
                Returns names of contained (geneModel) list.}
        %
        \item{names<-}{\code{signature(object= "geneList",i="numeric"}
                Sets names of contained (geneModel) list.}
        %
        \item{+}{\code{signature(e1= "geneList", e2="geneList"}
            Combines two lists together.}
        %
    }
}
%
\usage{
geneList(ref, genes, interior=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{ref}{ensemblGenome or ucscGenome.
        Genome object from which gene data is extracted.}
    %
    \item{genes}{Vector of gene_id's}
    %
    \item{interior}{When \code{FALSE}, exon and transcript data
        skipped.}
}
%
\author{
Wolfgang Kaisers
}
\examples{
##-------------------------------------##
## A) Ensembl
##-------------------------------------##
ensfile <- system.file("extdata",
                        "hs.ensembl.62.small.RData", package="refGenome")
ens <-loadGenome(ensfile)
gt <- getGeneTable(ens)
gl <- geneList(ens, gt$gene_id, interior=TRUE)
names(gl)
length(gl)
gl
gl[1]
gl[1] + gl[2]
}
\keyword{classes}
