\name{overlap}
\alias{overlap}
\title{ overlap function }
\description{
Overlaps query ranges with reference ranges. The function assumes that there is no overlap between reference ranges.
}
\usage{
overlap(qry, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qry}{data.frame with query ranges. qry should have columns 'id','start' and 'end'. The routine assumes that the table is ascending sorted by column 'start'.}
  \item{ref}{data.frame with reference ranges. ref should have columns 'id','start' and 'end'. The routine assumes that the table is ascending sorted by column 'start'.}
}
\details{
The routine assumes that qry and ref tables are ascending sorted by column 'start'. Otherwise the result will be incorrect (i.e. missing hits). The function assumes that there is no overlap between reference ranges. It will otherwise return only one, possibly arbitrary, hit per query range.
}
\value{
The function returns a data.frame
  \item{overlap}{Factor which encodes type of overlap between query and reference range. Levels no (no overlap; refid is set to 0),l (qry left overhangs ref), n (qry is contained in ref ), r (qry right overhangs ref)}
  \item{leftDiff}{Distance on left side between margins of query and reference.}
  \item{rightDiff}{Distance on right side between margins of query and reference.}
  \item{queryid}{id from qry table.}
  \item{refid}{id from ref table.}
}
\author{
Wolfgang Kaisers
}
\examples{
qry<-data.frame(id=1:6,start=c(10,18,61,78,82,110),end=c(15,22,63,87,90,120))
ref<-data.frame(id=1:5,start=c(20,40,60,80,100),end=c(25,45,65,85,105))
overlap(qry,ref)
}

\keyword{overlap}
