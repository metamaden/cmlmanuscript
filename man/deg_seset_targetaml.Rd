\name{deg_seset_targetaml}
\alias{deg_seset_targetaml}
\docType{data}
\title{
TARGET pediatric AML SummarizedExperiment object
}
\description{
Annotated gene expression data for TARGET pediatric AML dataset, in SummarizedExperiment format.
}
\usage{data("deg_seset_targetaml")}
\format{
  A SummarizedExperiment object of gene expression data for the TARGET pediatric AML dataset.
}
\details{
}
\source{
}
\references{
(website) TARGET. <https://ocg.cancer.gov/programs/target>
}
\examples{
data(deg_seset_targetaml)

dim(deg.seset)
# [1] 1984  145

class(deg.seset)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

dim(colData(deg.seset))
# [1] 145  85

}
\keyword{datasets}
\keyword{gene expression}
\keyword{pediatric AML}
\keyword{SummarizedExperiment}
