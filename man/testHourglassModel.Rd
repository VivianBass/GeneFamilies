% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression_funks.R
\name{testHourglassModel}
\alias{testHourglassModel}
\title{Test the hypothesis of the developmental hourglass by carrying out various
t-tests. In these, compare the distributions of expression profile distances
within a given tissue between gene classes. This should show which gene
class contributes to the hourglass model. Also compare for a given gene
class early and late developmental tissues with intermediate samples. This
should show wether early and late stages are more diverse than intermediate
ones, and thus give the hourglass shape.}
\usage{
testHourglassModel(x, tissues = c("seedling", "cotyledon",
  "developing leaf", "flower stage 9", "flower stage 16"),
  gene.classes = c("Orthologs", "Non-Orthologs"),
  early.late.indx = c(2, 4), intermediate.indx = 3,
  statistic = "median")
}
\arguments{
\item{x}{an instance of \code{base::data.frame} with columns: 'Family',
'Statistic', 'Gene.Class', 'Value', 'Tissue'. For example this package's
data \code{tandems.exp.prof.dists.tissue.orth.dist}.}

\item{tissues}{a character vector identifying the tissues to be
investigated. Default is \code{c('seedling', 'cotyledon', 'developing leaf',
'flower stage 9', 'flower stage 16')}.}

\item{gene.classes}{a character vector identifying the gene classes to be
investigated. Default is \code{c('Orthologs', 'Non-Orthologs')}.}

\item{early.late.indx}{A numeric vector identifying which entries of
\code{tissues} are considered early and late developmental stages. Default
is \code{c(2,4)}.}

\item{intermediate.indx}{A single integer identifying which entry in
\code{tissues} is considered to be the intermediate developmental stage.
Default is \code{3}.}

\item{statistic}{a string identifying which statistic of \code{x} to be
investigated. Default is \code{'median'}.}
}
\value{
A \code{base::data.frame} with columns: 'Tissue',
'Alternative.Hypothesis', 'p.value', 'effect.size', 'p.value.adjusted'
}
\description{
Test the hypothesis of the developmental hourglass by carrying out various
t-tests. In these, compare the distributions of expression profile distances
within a given tissue between gene classes. This should show which gene
class contributes to the hourglass model. Also compare for a given gene
class early and late developmental tissues with intermediate samples. This
should show wether early and late stages are more diverse than intermediate
ones, and thus give the hourglass shape.
}
