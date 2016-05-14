\name{transcriptModel-class}
\Rdversion{1.1}
\docType{class}
\alias{transcriptModel-class}
\alias{initialize,transcriptModel-method}
%
\alias{getTranscript}
\alias{getTranscript-methods}
\alias{getTranscript,geneModel,numeric-method}
\alias{getTranscript,geneModel,character-method}
%
\alias{getExonData}
\alias{getExonData-methods}
\alias{getExonData,transcriptModel-method}
%
\alias{getCdsData}
\alias{getCdsData-methods}
\alias{getCdsData,transcriptModel-method}
%
\alias{plot.transcriptModel}
%
\title{Class \code{"transcriptModel"}}
\description{
transcriptModel represents data for a single gene.

}
%
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{.transcriptModel()}.
}
%
\section{Slots}{
    \describe{
        \item{\code{id}:}{Object of class \code{"character"}.
            ID of one single transcript.}
        %
        \item{\code{name}:}{Object of class \code{"character"}.
            Name of one single transcript.}
        %
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
        \item{\code{biotype}:}{Object of class \code{"character"}.
            Biotype of transcript.}
        %
        \item{\code{coords}:}{Object of class \code{"character"}.
            Limiting coordinates (start, end) of gene.}
        %
        \item{\code{exons}:}{Object of class \code{"data.frame"}
            with two columns.
            Contains limiting coordinates (start, end) of exons.}
        %
        \item{\code{cds}:}{Object of class \code{"data.frame"}
            with two columns.
            Contains limiting coordinates (start, end) of CDS.}
        %
        \item{\code{stcodon}:}{Object of class \code{"integer"}.
            Coordinates start and stop codons.}
        %
        \item{\code{prime_utr}:}{Object of class \code{"integer"}.
            Coordinates five- and three-prime_utr.}
        %
        \item{\code{version}:}{Object of class \code{"integer"}.
            Version as stated in ensembl GTF files.}
    }
}

\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "transcriptModel")}:
            Creates a sensible printout. }
    }
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
gm <- geneModel(ens, gene_id)
tr <- getTranscript(gm, 1)
tr
getExonData(tr)
getCdsData(tr)
# + + + + + + + + + + + + + + + + + + #
# B) UCSC
# + + + + + + + + + + + + + + + + + + #
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
gt <- getGeneTable(uc)
gene_id <- as.character(gt$gene_id[1])
gm <- geneModel(uc, gene_id)
tr <- getTranscript(gm, 1)
tr

}
\keyword{classes}
