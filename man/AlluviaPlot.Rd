% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{AlluviaPlot}
\alias{AlluviaPlot}
\title{Plot Alluvial Plot}
\usage{
AlluviaPlot(
  cellmeta,
  by,
  fill,
  colors = NULL,
  bar.width = 0.5,
  legend.ncol = 1
)
}
\arguments{
\item{cellmeta}{a data frame containing the cell metadata.}

\item{by}{a character string specifying the column name in \code{cellmeta} which will be used to group the cells.}

\item{fill}{a character string specifying the column name in \code{cellmeta} which will be used to fill the plot.}

\item{colors}{(optional) vector of colours to use for filling the plot. If not specified, the default colour scheme will be used.}

\item{bar.width}{a numeric value between 0 and 1 specifying the width of the bars in the plot.}

\item{legend.ncol}{an integer specifying the number of columns in the legend.}
}
\value{
A ggplot object representing an alluvial plot.
}
\description{
Plot Alluvial Plot
}
\examples{
AlluviaPlot(mtcars, "cyl", "gear")

}
\references{
https://stackoverflow.com/questions/73372641/shaded-area-between-bars-using-ggplot2
}
\concept{plot}
