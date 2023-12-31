% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/BN-GWAS_applications_get_gene_gene_network_from_GTEx.R
\name{get_Residualized_GTEx}
\alias{get_Residualized_GTEx}
\title{Prepare EXP Matrix}
\usage{
get_Residualized_GTEx(
  TPM_thres,
  low_var_thres,
  tissue_type_to_study,
  attribute,
  finaldf,
  pca_res,
  subject_pheno
)
}
\arguments{
\item{TPM_thres}{parameter used to filter genes with lower TPM}

\item{low_var_thres}{parameter used to filter genes with low variance}

\item{tissue_type_to_study}{target tissue from GTEx}

\item{attribute}{attributes file from GTEx recording all related attributes of subjects}

\item{finaldf}{GTEx RNA-seq data}

\item{pca_res}{PCA analysis results for all recruited subjects}

\item{subject_pheno}{subject phenotype data}
}
\value{
resid_mat returned adjusted tissue-specific RNA-seq data
}
\description{
The aim of this function is to prepare the gene expression matrix, controlling for possible batch effects (and other confounding factors)
}
\examples{
get_Residualized_GTEx(TPM_thres, low_var_thres, tissue_type_to_study, attribute, finaldf, pca_res, subject_pheno)

}
