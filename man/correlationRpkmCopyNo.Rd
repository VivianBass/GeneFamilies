% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression_funks.R
\name{correlationRpkmCopyNo}
\alias{correlationRpkmCopyNo}
\title{Generates a data.frame of the Genes' names, copy numbers, and expression
values (RPKM). Furthermore, computes a lm model for the observed data. If
expression levels are extracted for more than a single tissue the
\code{mean} value is used. Note, that only genes are considered that have
some expression in any tissue, never expressed genes are discarded from the
analysis.}
\usage{
correlationRpkmCopyNo(rpkm.df = rpkm.rna.seq.counts, r.gene.col = "id",
  r.rpkm.col = "expression", r.tissue.col = "tissue",
  copy.no.df = gene.copy.number.df, c.gene.col = "Gene.no.expr.var",
  c.copy.no.col = "copy.no", tissues = unique(rpkm.df[, r.tissue.col]),
  stat.funk = mean, lapply.funk = lapply)
}
\arguments{
\item{rpkm.df}{An instance of data frame holding the genes expression counts
(RPKM). Default is \code{rpkm.rna.seq.counts}.}

\item{r.gene.col}{A string or integer identifying the column of
\code{rpkm.df} in which to find the gene accessions (IDs). Default is
\code{'id'}.}

\item{r.rpkm.col}{A string or integer identifying the column of
\code{rpkm.df} in which to find the genes' expression counts (RPKM). Default
is \code{'expression'}.}

\item{r.tissue.col}{A string or integer identifying the column of
\code{rpkm.df} in which to find the tissue in which the expression was
measured. Default is \code{'tissue'}.}

\item{copy.no.df}{An instance of \code{base:data.frame} in which to find the
genes' copy numbers. Default is \code{'gene.copy.number.df'}.}

\item{c.gene.col}{A string or integer identifying the column of
\code{copy.no.df} in which the gene accessions (IDs) are stored. Default is
\code{'Gene.no.expr.var'}.}

\item{c.copy.no.col}{A string or integer identifying the column of
\code{copy.no.df} in which the genes' copy numbers are held. Default is
\code{'copy.no'}.}

\item{tissues}{A character of minimum length one in which the tissues to
extract expression levels for are stored. Default is \code{unique(rpkm.df[,
r.tissue.col])}}

\item{stat.funk}{The statistical measure applied on the expression values
found for a respective gene within the argument \code{tissues}. Default is
\code{base::mean}.}

\item{lapply.funk}{One of \code{base::lapply} or \code{parallel::mclapply}
to indicate which function shall be used to iterate over the respective
genes. Default is \code{lapply}.}
}
\value{
A list with entries: \code{data} the above data.frame, \code{lm} the
above mentioned computed generalized linear model, and the \code{R^2} (R
squared) value.
}
\description{
Generates a data.frame of the Genes' names, copy numbers, and expression
values (RPKM). Furthermore, computes a lm model for the observed data. If
expression levels are extracted for more than a single tissue the
\code{mean} value is used. Note, that only genes are considered that have
some expression in any tissue, never expressed genes are discarded from the
analysis.
}
