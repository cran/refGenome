\name{refGenome-class}
\Rdversion{1.1}
\docType{class}
\alias{refGenome-class}
\alias{initialize,refGenome-method}
\alias{basedir}
\alias{basedir-methods}
\alias{basedir,refGenome-method}
\alias{basedir<-}
\alias{basedir<--methods}
\alias{basedir<-,refGenome-method}
\alias{extractFeature}
\alias{extractFeature-methods}
\alias{extractFeature,refGenome-method}
\alias{extractSeqids}
\alias{extractSeqids-methods}
\alias{extractSeqids,refGenome-method}
\alias{extractTranscript}
\alias{extractTranscript-methods}
\alias{extractTranscript,refGenome-method}
\alias{getGtf}
\alias{getGtf-methods}
\alias{getGtf,refGenome-method}
\alias{setGtf}
\alias{setGtf-methods}
\alias{setGtf,refGenome-method}
\alias{getAttr}
\alias{getAttr-methods}
\alias{getAttr,refGenome-method}
\alias{setAttr}
\alias{setAttr-methods}
\alias{setAttr,refGenome-method}
\alias{read.gtf}
\alias{read.gtf-methods}
\alias{read.gtf,refGenome-method}
\alias{tableAttributeTypes}
\alias{tableAttributeTypes-methods}
\alias{tableAttributeTypes,refGenome-method}
\alias{tableSeqids}
\alias{tableSeqids-methods}
\alias{tableSeqids,refGenome-method}
\alias{ucPrimAssembly}
\alias{ensPrimAssembly}

\title{Class \code{"refGenome"}}
\description{
refGenome class: Virtual base class for ucscGenome and ensemblGenome.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{refGenome(dbfile)}.
dbfile represents SQLite database file.
}
\section{Slots}{
    \describe{
        \item{\code{basedir}:}{Object of class \code{"character"}
            Directory where SQLite database is written.}
        %
        \item{\code{ev}:}{Object of class \code{"environment"}
            Environment that contains data structures.
            Optionally, there are gtf and attr data.frames. }
        %
  }
}
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "refGenome")}:
            Creates a sensible printout.}
        %
        \item{tableSeqids}{\code{signature(object = "refGenome")}:
            Returns tabled seqids which counts all annotations for each seqid.
            The regex argument will display, which seqid is covered by regex. 
            Intended as preparation for extractSeqids.}
        %
        \item{extractSeqids}{\code{signature(object = "refGenome")}:
            Returns a filtered version of object.
            Only data for the 'regex' specified seqids is contained. 
            The functions 'ucPrimAssembly' and 'ensPrimAssembly' return regular
            expressions that allow extraction of primary assemblies for UCSC
            and Ensembl respectively.}
        %
        \item{getGtf}{\code{signature(object = "refGenome")}:
            Returns content of gtf table.}
        %
        \item{setGtf}{\code{signature(object = "refGenome")}:
            Writes content of gtf table.}
        %
        \item{getAttr}{\code{signature(object = "refGenome")}:
            Returns content of attribute table.}
        %
        \item{setAttr}{\code{signature(object = "refGenome")}:
            Writes content of attribute table.}
        %
        \item{read.gtf}{\code{signature(object = "refGenome")}: 
            Imports content of gtf file. 
            This is the basic mechanism for data import. 
            It works the same way for ucscGenome and for ensemblGenome. 
            The attribute items are parsed and written in parts to gtf table
            and attribute table.}
        %
        \item{tableAttributeTypes}{\code{signature(object= "refGenome")}:
            Prints out number of rows for each type in attribute table.}
        %
        \item{saveGenome}{\code{signature(object = "refGenome",
            filename="character", useBasedir="logical")}: 
            Saves content of all tables to RData file. 
            When useBasedir is set (default), basedir (from basedir-slot)
            is prefixed.}
        %
        \item{writeDB}{\code{signature(object = "refGenome")}:
            Copies content of gtf, attr and xref table to database. }
        %
    }
}
\author{Wolfgang Kaisers}
\examples{
# + + + + + + + + + + + + + + + + + + #
# Loading sample data:
# + + + + + + + + + + + + + + + + + + #
ensfile<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)
ens
ddx<-extractByGeneName(ens,"DDX11L1")
ddx
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
uc
ddx<-extractByGeneName(uc,"DDX11L1")
ddx
}
\keyword{classes}
