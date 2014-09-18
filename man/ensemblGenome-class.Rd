\name{ensemblGenome-class}
\Rdversion{1.1}
\docType{class}
\alias{ensemblGenome-class}
\alias{ensemblGenome}
%
\alias{extractPaGenes}
\alias{extractPaGenes-methods}
\alias{extractPaGenes,ensemblGenome-method}
%
\alias{moveAttributes}
\alias{moveAttributes-methods}
\alias{moveAttributes,ensemblGenome-method}
%
\alias{tableTranscript.id}
\alias{tableTranscript.id-methods}
\alias{tableTranscript.id,ensemblGenome-method}
%
\alias{tableFeatures}
\alias{tableFeatures-methods}
\alias{tableFeatures,ensemblGenome-method}
%
\alias{tableTranscript.name}
\alias{tableTranscript.name-methods}
\alias{tableTranscript.name,ensemblGenome-method}
%
\alias{getGeneTable}
\alias{getGeneTable-methods}
\alias{getGeneTable,ensemblGenome-method}
%

\title{Class \code{"ensemblGenome"}}
\description{
ensemblGenome represents ensembl genomic annotation data.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{ensemblGenome(dbfile)}.
    'dbfile' represents SQLite database file.
}
\section{Slots}{
    \describe{
        \item{\code{basedir}:}{Object of class \code{"character"} 
                        Directory where SQLite database is written.}
        %
        \item{\code{ev}:}{Object of class \code{"environment"} Environment that
        contains data structures. Optionally, there are gtf and attr
        data.frames. }
        %
    }
}
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "refGenome")}:
            Creates a sensible printout.}
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
        \item{getGeneTable}{\code{signature(object = "refGenome")}:
            Returns content of genes table when table exists. 
            Otherwise NULL is returned.}
        %
        \item{setAttr}{\code{signature(object = "refGenome")}: 
            Writes content of attribute table.}
        %
        \item{read.gtf}{\code{signature(object = "refGenome")}:
            Imports content of gtf file.
            This is the basic mechanism for data import. It works the
            same way for ucscGenome and for ensemblGenome. The attribute items
            are parsed and written in parts to gtf table and attribute table.}
        %
        \item{extractPaGenes}{\code{signature(object="ensemblGenome")}: 
            Extracts all annotations on primary assembly. The function returns
            a data.frame. Used as shortcut to directly extract a table from gtf
            files.}
        %
        \item{extractFeature}{\code{signature(object="ensemblGenome")}:
            Extracts annotated positions which are classified as given 'feature'
            argument. Returns an 'ensemblGenome' object.}
        %
        \item{extractByGeneName}{\code{signature(object="ensemblGenome",
            geneNames="character")}:
            Extracts ensemblGenome object which contains table subsets.
            When none of the \code{geneNames} matches, 
            the function returns \code{NULL}.}
        %
        \item{extractTranscript}{\code{signature(object="ensemblGenome",
            transcripts="character")}:
            Extracts ensemblGenome object which contains table subsets}
        %
        \item{getGenePositions}{\code{signature(object="ucscGenome",
            force="logical")}:
            Extracts table with position data for whole genes
            (smallest exon start position and largest exon end position. A copy
            of the table will be placed inside the internal environment. Upon
            subsequent call only a copy of the contained table is returned
            unless \code{force=TRUE} is given. Upon \code{force=TRUE} new gene 
            positions are calculated regardless of existing tables.)}
        %
        \item{tableTranscript.name}{\code{signature(object="ensemblGenome")}:
            Extracts table object which contains tabled 'transcript_name' column
            of gtf table}
        %
        \item{tableTranscript.id}{\code{signature(object="ensemblGenome")}:
            Extracts table object which contains tabled 'transcript_id' column
            of gtf table}
        %
        \item{writeDB}{\code{signature(object = "refGenome")}:
            Copies content of gtf, attr and xref table to database.}
}
}
\references{
http://www.ensembl.org/info/data/ftp/index.html
http://mblab.wustl.edu/GTF22.html#fields
}
\author{
Wolfgang Kaisers
}
\examples{
# + + + + + + + + + + + + + + + + + + + #
# Create an instance from scratch
# Real data:
# ftp://ftp.ensembl.org/pub/release-70/gtf/homo_sapiens/Homo_sapiens.GRCh37.70.gtf.gz
# + + + + + + + + + + + + + + + + + + #
ens<-ensemblGenome()
basedir(ens)<-system.file("extdata",package="refGenome")
ens_gtf<-"hs.ensembl.62.small.gtf"
read.gtf(ens,ens_gtf)
tableAttributeTypes(ens)
moveAttributes(ens,c("gene_name","transcript_name","exon_number"))
# Loading of a saved genome:
ensfile<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)

# + + + + + + + + + + + + + + + + + + + # 
# Saving and loading 
# Save as R-image (fast loading)
# + + + + + + + + + + + + + + + + + + #
basedir(ens)<-getwd()
saveGenome(ens,"hs.ensembl.62.small.RData",useBasedir=FALSE)
enr<-loadGenome("hs.ensembl.62.small.RData")

# Save as SQLite database
# + + + + + + + + + + + + + + + + + + + #
# Commented out because RSQlite 
# seems to produce memory leaks
# + + + + + + + + + + + + + + + + + + + #
#writeDB(ens,filename="ens62.db3",useBasedir=FALSE)
#edb<-loadGenomeDb(filename="ens62.db3")

# + + + + + + + + + + + + + + + + + + + # 
# Extract data for Primary Assembly seqids
# + + + + + + + + + + + + + + + + + + #
enpa<-extractSeqids(ens,ensPrimAssembly())
# Tables all features in 'gtf' table
tableFeatures(enpa)
# Extract Coding sequences for Primary Assemblys
enpafeat<-extractFeature(enpa,"exon")
# Shortcut. Returns data.frame
engen<-extractPaGenes(ens)

# + + + + + + + + + + + + + + + + + + + # 
# Extract data for indival Genes
# + + + + + + + + + + + + + + + + + + #
ddx<-extractByGeneName(ens,"DDX11L1")
ddx
tableTranscript.id(ddx)
tableTranscript.name(ddx)
fam<-extractTranscript(ens,"ENST00000417324")
fam
# Extract range limits of entire Genes
gp<-getGenePositions(ens)
gp
# + + + + + + + + + + + + + + + + + + + #
}
\keyword{classes}
