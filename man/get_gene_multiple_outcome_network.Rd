% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/BN-GWAS_application_get_gene_outcome_causal_network.R
\name{get_gene_multiple_outcome_network}
\alias{get_gene_multiple_outcome_network}
\title{Infer Multiple Gene Outcome Causal Network}
\usage{
get_gene_multiple_outcome_network(
  outcome_name,
  expr,
  pheno,
  predicted_pheno,
  PCs,
  alpha,
  num_workers
)
}
\arguments{
\item{outcome_name}{updated coded ID for interested outcome variable(apply the inner_join function, the column name of outcome will change,e.g.,30690-0.0->X30690.0.0}

\item{expr}{genotype-predicted tissue-specific expression profiles from PrediXcan with subjects IDs in the 1st column named as FID}

\item{pheno}{target phenotype data extracted with subjects IDs in the 1st column named as FID}

\item{predicted_pheno}{PRS predicted phenotype(s) with subjects IDs in the 1st column named as FID, the column names of the predicted phenotype(s) should also start with "ENSG"}

\item{PCs}{Top 10 PCs with subjects IDs in the 1st column named as FID}

\item{alpha}{p-value threshold used to define causal relations}

\item{num_workers}{number of cores used for parallel computing}
}
\value{
gene-outcome causal network object
}
\description{
The aim of this function is to infer the gene outcome causal network for multiple outcomes.
}
\examples{
get_gene_multiple_outcome_network(outcome_name, expr, pheno, predicted_pheno, PCs, alpha,num_workers)

}
