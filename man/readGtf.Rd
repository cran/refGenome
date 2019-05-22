\name{read.gtf}
\Rdversion{1.1}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{read.gtf}
\alias{read.gtf-methods}
\alias{read.gtf,refGenome-method}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Title
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\title{Reading and parsing GTF files into refGenome objects.}
%
\description{Reads and parses content of GTF files. The file content is
written into the provided object into the environment located in 'ev' slot
(i.e. per reference). The function writes two tables: 'gtf' containing the main
file content and 'genes' containing data from 'gene' typed features.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\usage{
read.gtf(object, filename="transcripts.gtf", sep = "\t",
            useBasedir=TRUE, comment.char = "#", progress=100000L, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{\code{refGenome} object. Will contain the extracted
        data.}
    %
    \item{filename}{(Base-)Name of GTF file.}
    %
    \item{sep}{Character: Column separator in GTF file. Standard value is '\\t'.}
    %
    \item{useBasedir}{Logical: Shall basedir (from refGenome object)
        be appended to filename?}
    %
    \item{comment.char}{Character: Lines beginning with this character
        will be skipped.}
    %
    \item{progress}{Integer: The parsing routine prints a progress
        Information after reading the given number of lines.}
    %
    \item{...}{Currently unused.}
}
\details{GTF is an extension of the GFF file format. GTF contains tabled data:
Nine columns separated by a tab delimiter. The last column expands into a
list of attributes, separated by a semicolon an exactly one space. Each
attribute consists of a type - value pair which are separated by one empty
space. Enclosing quotation marks (") around attribute values are  marks are
skipped during import.
}
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
%% Returned value
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
\value{None. The provided object is filled with the parsed data. Two tables
are generated: 'gtf' and 'genes'. The first eight columns of the gtf table
are fixed. The content is described in the following table.
\tabular{ll}{
  id \tab Numeric index for unique site. Integer.\cr
  seqid \tab Chromosome identifier. Character.\cr
  source \tab Program which generated data.\cr
  feature \tab Feature type (e.g. 'exon', 'CDS'). Character.\cr
  start \tab Start position of feature (1-based). Integer.\cr
  end \tab End position of feature (inclusive). Integer.\cr
  score \tab Value between 0 and 1000 ("." for no score). Character.\cr
  strand \tab '+', '-' or '.'. Character.\cr
  frame \tab 0-2 for coding exons. '.' otherwise. Character.\cr
}
}
\author{
Wolfgang Kaisers
}
\references{UCSC Genome Bioinformatics: Data File Formats.
\url{http://genome.ucsc.edu/FAQ/FAQformat.html#format3}}
\examples{
##-------------------------------------##
## Ensembl
##-------------------------------------##
ef <- system.file("extdata", package="refGenome")
en <- ensemblGenome(ef)
read.gtf(en, "hs.ensembl.76.small.gtf")
}
\keyword{refJunctions}
