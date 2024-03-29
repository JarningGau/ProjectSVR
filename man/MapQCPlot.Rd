% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{MapQCPlot}
\alias{MapQCPlot}
\title{Visualizing the mapping quality of projected query cells.}
\usage{
MapQCPlot(seu.q, p.adj.cutoff = 1e-04, map.q.cutoff = NULL)
}
\arguments{
\item{seu.q}{The Seurat object, output of ProjectSVR::MapQuery().}

\item{p.adj.cutoff}{The adjusted p value cutoff. Default: 1e-4}

\item{map.q.cutoff}{The cutoff for mean.knn.dist. If not provided, the maximum value where mapQ.p.adj > 0 will be used.}
}
\value{
A combined plot for mapping quality.
}
\description{
Visualizing the mapping quality of projected query cells.
}
\concept{plot}
