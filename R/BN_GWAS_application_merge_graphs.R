#' Merge graph
#' 
#' Merge gene-gene and gene-outcome graph into a whole graph
#' 
#' @param gene_graph_adj m x m adjacency matrix for the gene-gene causal network: [i,j] =1 indicates an edge from i to j (graph package documentation)
#' @param gene_to_outcome m x 1 causal relations vector between genes and outcome: gene_to_outcome[i]=1 indicates causal connection from gene i to outcome
#'  
#' @return (m+1) x (m+1) transposed adjacency matrix bgadd2[i,j]=1 indicates an edge from j to i in the causal graph
#' 
#' @examples 
#' merge_causal_graphs(gene_graph_adj, gene_to_outcome)
#' 
#' @export 
merge_causal_graphs <- function(gene_graph_adj, gene_to_outcome){

  # *************************************
  # This function adds background knowledge to a graph from pcfit
  # ****************************************
  # library(graph)
  # library(pcalg)
  
  # *******************************************************************
  # get adjacency matrix for GTEx network derived from bootstrapped pc
  # *******************************************************************
  adj_mat <- gene_graph_adj
  final_mat <- cbind(adj_mat, gene_to_outcome) ## [i,j] means an edge from i to j (here it is gene -> disease, but NOT the other way round
  final_mat <- rbind(final_mat, rep(0, ncol(adj_mat) + 1)) ## need to check if the diagonals are all zero 
 
  # **********************************************************
  # Merge two graph by add background knowledge
  # **********************************************************
  # the gInput takes [i,j] to mean an edge from j to i, which is opposite to the convention in graph package
  mat_input <- t(final_mat) 
  t1 <- proc.time()
  bgadd2 <- addBgKnowledge(gInput = mat_input,  
                           verbose = TRUE, 
                           checkInput = FALSE)					   
  proc.time() - t1
  return(bgadd2)
}

#' Get IDA and joinytIDA
#' 
#' Estimate total and direct causal effects of genes on the outcome
#' 
#' @param estimateGraph estimated adacency matrix for the whole causal graph including both "genes" and phenotype(s)
#' @param cov_mat the covariance matrix for the genes and phenotypes
#' @param target_genes target genes for causal effects estimation
#'  
#' @return Effects_mat_final returned causal effects matrix with estimated causal effect 
#' 
#' @examples 
#' Get_IDA_and_jointIDA(estimateGraph, cov_mat, target_genes)
#' 
#' @export 
Get_IDA_and_jointIDA <- function(estimateGraph,cov_mat,target_genes){
  graph.fit <- as(estimateMatrix, "graphNEL") # Convert the combined graph as a grahNEL object
  nodes_estimateGraph <- nodes(estimateGraph)
  phenotype_loc <- length(nodes_estimateGraph)  # phenotype is the last node in the Graph
  x.pos <- nodes_estimateGraph[!nodes_estimateGraph %in% phenotype_loc]
  y.pos <- phenotype_loc
  x.pos <- as.numeric(x.pos) # change the data type from character to numeric
  y.pos <- as.numeric(y.pos)
  estimate_jointIDA_list <- jointIda(x.pos=x.pos,
                                     y.pos=y.pos,
                                     cov_mat, cov_mat,
                                     graphEst=estimateGraph, technique="RRC")
  cat("The estimation for jointIDA is completed!\n")
  jointIDAs_len <- ncol(estimate_jointIDA_list)
  jointIDAs_suffix <- seq(1,jointIDAs_len)
  jointIDAs_names <- paste0("Direct_effect_",jointIDAs_suffix)
  colnames(estimate_jointIDA_list) <- jointIDAs_names
  # Create a matrix to store the IDA and jointIDA for each causal genes
  Effects_mat <- matrix(0,nrow=length(x.pos),ncol=2)
  colnames(Effects_mat) <- c("Variable","Total_effect")
  
  for (i in 1: length(x.pos)){
    cat("This is the estimation of IDA for the",i,"variable.\n")
    node1 <- nodes_estimateGraph[i] # node1 is gene
    node2 <- phenotype_loc # node2 is phenotype in this case
    IDA.estimate <- ida(node1,node2, cov_mat, graph.fit, method = "local", type = "pdag")	
    Effects_mat[i,1] <- target_genes[i]
    # Revised by: yinly, June 23, 2020
    # Description: the length of the derived IDA.estimate may be larger than 1
    if(length(IDA.estimate) == 1){
      Effects_mat[i,2] <- IDA.estimate
    } else {
      IDA.estimate_combine <- NULL
      for (j in 1: length(IDA.estimate)){
        if (is.null(IDA.estimate_combine)){
          IDA.estimate_combine <- IDA.estimate[j]
        } else {
          IDA.estimate_combine <- paste(IDA.estimate_combine,IDA.estimate[j],sep = ",")
        }
      }
      Effects_mat[i,2] <- IDA.estimate_combine
    }   
  }
  cat("The estimation for IDA is completed!\n")
  Effects_mat_final <- cbind(Effects_mat,estimate_jointIDA_list)
  return(Effects_mat_final)
}

#' Pvalue-adjusted Covariance Matrix  
#' 
#' Get pvalue-adjusted covariance matrix
#' 
#' @param rho_pMax m x m adjacency matrix for the gene-gene causal network: [i,j] =1 indicates an edge from i to j (graph package documentation)
#' @param cov calcuated covariance matrix between genes and outcome(s)
#'  
#' @return returned adjusted covariance matrix
#' 
#' @examples 
#' get_pvalue_adjusted_cov(rho_pMax, cov)
#' 
#' @export 
get_pvalue_adjusted_cov <- function(rho_pMax, cov){
  diag(rho_pMax) <- 0
  rho_pMax[rho_pMax == -1] <- 1
  cov_adjusted <- glassoFast(cov, rho = rho_pMax)$w
  return(cov_adjusted)
}

#' Iterative Pvalue-adjusted Covariance Matrix  
#' 
#' get iterative pvalue-adjusted covariance matrix  
#' 
#' @param rho_pMax m x m adjacency matrix for the gene-gene causal network: [i,j] =1 indicates an edge from i to j (graph package documentation)
#' @param cov calcuated covariance matrix between genes and outcome(s)
#' @param iterateNum number of iterations
#' @param thres threshold used to determine whether the penalization process is coverged!
#' @param n sample size of the pheno dataset used for covariance matirx calculation
#'  
#' @return returned adjusted covariance matrix
#' 
#' @examples 
#' get_iterative_pvalue_adjusted_cov(rho_pMax, cov, iterateNum, thres)
#' 
#' @export 
get_iterative_pvalue_adjusted_cov <- function(rho_pMax, cov, iterateNum, thres, n){
  diag(rho_pMax) <- 0
  rho_pMax[rho_pMax==-1] <- 1
  iteration_error <- array(0, iterateNum)
  thresh <- thres
  cov_pvalue_iterative_adjusted <- cov
  for (i in 1:iterateNum){
    GLASSO <- glassoFast(cov_pvalue_iterative_adjusted, rho=rho_pMax)
    cov_pvalue_iterative_adjusted <- GLASSO$w
    iteration_error[i] <- iteration_error[i] + sum(GLASSO$wi[upper.tri(GLASSO$wi)] != 0) *
      log(n)/2
    if ((i>1) && (abs(iteration_error[i]-iteration_error[i-1])< thresh)){
      break
    }
  }
  return(cov_pvalue_iterative_adjusted)
}