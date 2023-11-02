## -----------------------------------------------------------------------------
set.seed(123)

library(BNGWAS)
library(data.table)

## -----------------------------------------------------------------------------
resid_mat <- fread("./example_resid_mat.csv")
n.cores <- 15
glasso_obj <- get_Glasso_Graph(resid_mat,n.cores)

alpha <- 0.05
numCores <- 15
m.max <- 3
pc_obj <- get_Bayesian_Network(alpha,numCores,m.max,resid_mat,glasso_obj)
print("Finished studying the gene-gene causal network")

## -----------------------------------------------------------------------------
expr <- fread("./demo_expr.csv")
pheno<- fread("./demo_pheno.csv")
PCs<- fread("./demo_PCs.csv")
predicted_pheno <- fread("./demo_predicted_pheno.csv")
num_workers <- 10
outcome_name<-"BMI"
pcSimple.fit <- get_gene_single_outcome_network(outcome_name,expr,pheno,PCs,alpha,num_workers)
pcSimple.fit.multiple <- get_gene_multiple_outcome_network(outcome_name,expr,pheno,predicted_pheno,PCs,alpha,num_workers)
print("Finished studying the gene-outcome causal network")

## -----------------------------------------------------------------------------
pc_graph_adj <- pc_obj@graph
gene_graph_adj <- as(pc_graph_adj, "matrix")

genes_names <- names(pcSimple.fit$G)
assoc <- as.vector(unlist(pcSimple.fit$G))
zMin <- as.vector(unlist(pcSimple.fit$zMin))
gene_outcome_network <- as.data.frame(cbind(genes_names,assoc,zMin))
gene_to_outcome <-rep(0,length(assoc))
assoc_index <- which(gene_outcome_network$assoc==TRUE)
gene_to_outcome[assoc_index] <- 1
bgadd2 <- merge_causal_graphs(gene_graph_adj,gene_to_outcome)

outcome <- pheno
genes_trait <- inner_join(expr,outcome, by="FID")
gene_trait_valid <- genes_trait[,-1]
cat("The sample size is",nrow(gene_trait_valid),"\n")
cov_mat <- cov(gene_trait_valid)
target_genes <- colnames(expr)[-1]

pheno_num <- 1
gene_num <- length(gene_to_outcome)
rho_pMax_len <- gene_num + pheno_num
rho_pMax <- matrix(0,nrow= rho_pMax_len,ncol=rho_pMax_len)

# Get rho_pMax matrix for the overall causal graph
rho_pMax_gene <- pc_obj@pMax ##extract pMax matrix for gene-gene causal network
graph_estimate <- pc_obj@graph ## extract the inferred gene-gene causal graph
adj_temp <- as(graph_estimate,"matrix") ##convert the inferred gene-gene causal graph to adjacency matrix
rho_pMax_gene[which(rho_pMax_gene<alpha & adj_temp ==0)] = 1  

rho_pMax_gene_outcome <- 2 * pnorm(-abs(as.numeric(as.vector(unlist(gene_outcome_network$zMin))))) ##convert zMin for gene-outcome causal network to pvalues
rho_pMax[1:gene_num,1:gene_num] <- rho_pMax_gene # 
rho_pMax[,rho_pMax_len] <- c(rho_pMax_gene_outcome,-1)
rho_pMax[rho_pMax_len,] <- c(rho_pMax_gene_outcome,-1)
diag(rho_pMax) <- 0
rho_pMax[rho_pMax==-1] <- 1

estimateMatrix <- bgadd2  
estimateMatrix <- t(estimateMatrix) ## to make bgadd2 consistent with the graph derived from GTEx
rownames(estimateMatrix) <- seq(1,rho_pMax_len)
colnames(estimateMatrix) <- seq(1,rho_pMax_len)
# Convert the combined graph as a grahNEL object
graph.fit <- as(estimateMatrix, "graphNEL")
#convert format of the graph that have negative edges
estimateGraph <- graph.adjacency(estimateMatrix, mode="directed", weighted=TRUE)
estimateGraph <- igraph.to.graphNEL(estimateGraph)
Effects_mat_final <- Get_IDA_and_jointIDA(estimateGraph,cov_mat,target_genes)

cov_adjusted <- get_pvalue_adjusted_cov(rho_pMax,cov_mat)
Effects_mat_pval_adjusted_final <- Get_IDA_and_jointIDA(estimateGraph,cov_adjusted,target_genes)

iterateNum <-20 # maximum iteration number
thres <- 1e-04 # threshold used to determine whether the penalizaiton converged
n <- nrow(gene_trait_valid)
cov_iterative_adjust <- get_iterative_pvalue_adjusted_cov(rho_pMax,cov_mat,iterateNum,thres,n)
print("Finished merging the gene-gene and gene-outcome causal network")

