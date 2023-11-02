#' Infer Single Gene Outcome Causal Network
#' 
#' The aim of this function is to infer the gene outcome causal network for a single outcome.
#' 
#' @param outcome_name updated coded ID for interested outcome variable(apply the inner_join function, the column name of outcome will change,e.g.,30690-0.0->X30690.0.0
#' @param expr genotype-predicted tissue-specific expression profiles from PrediXcan with subjects IDs in the 1st column named as FID
#' @param pheno target phenotype data extracted with subjects IDs in the 1st column named as FID
#' @param PCs Top 10 PCs with subjects IDs in the 1st column named as FID
#' @param alpha p-value threshold used to define causal relations
#' @param num_workers number of cores used for parallel computing
#'  
#' @return gene-outcome causal network object
#' 
#' @examples 
#' get_gene_single_outcome_network(outcome_name, expr, pheno, PCs, alpha, num_workers)
#' 
#' @export 
get_gene_single_outcome_network <- function(outcome_name, expr, pheno, PCs, alpha, num_workers){
    outcome <- data.frame(pheno)

    # *************************************
    # merging expression with the outcome
    # *************************************
    # Here, after apply the inner_join function, the column name of outcome changed to "X30690.0.0"(original:30690-0.0)
    df <- inner_join(outcome, expr, by="FID")   #merge(outcome, expr, by="FID")
    print((dim(df)))
    df <- inner_join(df,PCs,by ="FID")
    
    # *****************
    # count the no. of genes available
    # ******************
    grep_res <- grep("ENSG", colnames(df) )
    start.ind <- min(grep_res)
    end.ind <- max(grep_res)
    no_genes <- end.ind - start.ind + 1   
    
    # ********************************
    # Residualized X
    # *********************************
    resid.mat <- matrix(nrow = nrow(df), ncol= no_genes )
    attach(df)
    for (i in 1:no_genes) {
        lm.obj <- lm(df[,start.ind+ i-1] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
        resid.mat[,i] <- lm.obj$resid
    }
    colnames(resid.mat) <- colnames(df)[start.ind:end.ind]
    
    # ********************************
    # Residualized Y
    # *********************************
    # Revised by yinly on Jan 7,2020,the inner_join function changed the column name of the target_outcome_name
    lm.obj_y <- lm(df[,outcome_name] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
    resid_y <- lm.obj_y$resid
    detach(df)
    
    # Revised by yinly, June 4, 2021
    # Description: remove rankNorm, since this maybe change the distribution of the data
    dfnew <- data.frame(outcome= resid_y,  resid.mat)
    dfnew_rankNorm <- dfnew  
    # library(coop)
    expr_corMat <- pcor(dfnew_rankNorm[,-1])
    # now add back the correlation between the outcome and each feature
    # cf https://stackoverflow.com/questions/20410768/how-to-correlate-one-variable-to-all-other-variables-on-r
    cor_with_outcome <- apply(dfnew_rankNorm[,-1],  2 , function(col) {pcor(dfnew_rankNorm[,"outcome"], col)}  )
    precompute_corMat <- cbind(cor_with_outcome, expr_corMat)
    precompute_corMat <- rbind( c(1,cor_with_outcome), precompute_corMat ) 
    
    # ********************************************************
    # Apply PC-simple algorithm around the response variable on NFBC
    # **************************************************************
    # library(pcalg)   
    ## Load  data
    n <- nrow (dfnew_rankNorm)
    V <- colnames(dfnew_rankNorm) # labels aka node names  
    t1 <- proc.time()
    pcSimple.fit <- PCSelect_Parallel(y=dfnew_rankNorm[,1],
                                      dm=dfnew_rankNorm[,-1],
                                      method = c("parallel"),
                                      mem.efficient = FALSE,
                                      num_workers = num_workers,
                                      alpha=alpha,
                                      corMat = precompute_corMat,
                                      max_ord=3,
                                      corMethod = "standard",
                                      verbose = TRUE, directed = TRUE)
    proc.time()-t1
    return(pcSimple.fit)
}

#' Infer Multiple Gene Outcome Causal Network
#' 
#' The aim of this function is to infer the gene outcome causal network for multiple outcomes.
#' 
#' @param outcome_name updated coded ID for interested outcome variable(apply the inner_join function, the column name of outcome will change,e.g.,30690-0.0->X30690.0.0
#' @param expr genotype-predicted tissue-specific expression profiles from PrediXcan with subjects IDs in the 1st column named as FID
#' @param pheno target phenotype data extracted with subjects IDs in the 1st column named as FID
#' @param predicted_pheno PRS predicted phenotype(s) with subjects IDs in the 1st column named as FID, the column names of the predicted phenotype(s) should also start with "ENSG"
#' @param PCs Top 10 PCs with subjects IDs in the 1st column named as FID
#' @param alpha p-value threshold used to define causal relations
#' @param num_workers number of cores used for parallel computing
#'  
#' @return gene-outcome causal network object
#' 
#' @examples 
#' get_gene_multiple_outcome_network(outcome_name, expr, pheno, predicted_pheno, PCs, alpha,num_workers)
#' 
#' @export 
get_gene_multiple_outcome_network <- function(outcome_name, expr, pheno, predicted_pheno, PCs, alpha, num_workers){
    outcome <- data.frame(pheno)
    expr_new <- inner_join(expr, predicted_pheno,by="FID")
    
    # ************************************
    # merging expression with the outcome
    # ************************************
    expr <- data.frame(expr_new)
    # Here, after apply the inner_join function, the column name of outcome changed to "X30690.0.0"(original:30690-0.0)
    df <- inner_join(outcome, expr, by="FID")   #merge(outcome, expr, by="FID")
    print((dim(df)))
    df <- inner_join(df,PCs,by ="FID")

    # *****************
    # count the no. of genes available
    # ******************
    grep_res <- grep("ENSG", colnames(df))
    start.ind <- min(grep_res)
    end.ind <- max(grep_res)
    no_genes <- end.ind - start.ind + 1   
    
    # ********************************
    # Residualized X
    # *********************************
    resid.mat <- matrix(nrow = nrow(df), ncol= no_genes )
    attach(df)
    for (i in 1:no_genes) {
        lm.obj <- lm(df[,start.ind+ i-1] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
        resid.mat[,i] <- lm.obj$resid
    }
    colnames(resid.mat) <- colnames(df)[start.ind:end.ind]

    # ********************************
    # Residualized Y
    # *********************************
    # Revised by yinly on Jan 7,2020,the inner_join function changed the column name of the target_outcome_name
    lm.obj_y <- lm(df[,outcome_name] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
    resid_y <- lm.obj_y$resid
    detach(df)
    
    # Revised by yinly, June 4, 2021
    # Description: remove rankNorm, since this maybe change the distribution of the data
    dfnew <- data.frame(outcome= resid_y,  resid.mat)
    dfnew_rankNorm <- dfnew  
    # library(coop)
    expr_corMat <- pcor(dfnew_rankNorm[,-1])
    #now add back the correlation between the outcome and each feature
    #cf https://stackoverflow.com/questions/20410768/how-to-correlate-one-variable-to-all-other-variables-on-r
    cor_with_outcome <- apply(dfnew_rankNorm[,-1],  2 , function(col){pcor(dfnew_rankNorm[,"outcome"], col)})
    precompute_corMat <- cbind(cor_with_outcome, expr_corMat)
    precompute_corMat <- rbind(c(1,cor_with_outcome), precompute_corMat) 
    # ********************************************************
    # Apply PC-simple algorithm around the response variable on NFBC
    # **************************************************************
    # library(pcalg)   
    ## Load  data
    n <- nrow(dfnew_rankNorm)
    V <- colnames(dfnew_rankNorm) # labels aka node names  
    t1 <- proc.time()
    pcSimple.fit <- PCSelect_Parallel(y = dfnew_rankNorm[,1],
                                      dm = dfnew_rankNorm[,-1],
                                      method = c("parallel"),
                                      mem.efficient = FALSE,
                                      num_workers = num_workers,
                                      alpha = alpha,
                                      corMat = precompute_corMat,
                                      max_ord=3,
                                      corMethod = "standard",
                                      verbose = TRUE, directed = TRUE)
    proc.time()-t1
    return(pcSimple.fit)
}

#' 
#' 
#' Find columns in dm, that have nonzero parcor with y given any other set of columns in dm
#' 
#' @param y Response Vector (length(y)=nrow(dm))
#' @param dm Data matrix (rows: samples, cols: nodes)
#' @param method Character string specifying method; the default, "parallel" provides an parallelised method to implement all the conditional independence tests.
#' @param mem.efficient If TRUE, uses less amount of memory at any time point while running the algorithm
#' @param num_workers The numbers of cores CPU to run the algorithm
#' @param alpha Significance level of individual partial correlation tests
#' @param corMat Correlation matrix of the input variables
#' @param corMethod "standard" or "Qn" for standard or robust correlation estimation
#' @param max_ord Maximum number of variables for the conditional set
#' @param verbose 0-no output, 1-small output, 2-details
#' @param directed Logical; should the output graph be directed?
#'  
#' @return a list containing two object
#'      G: boolean vector with connected nodes and zMin: Minimal z values
#' 
#' @examples 
#' PCSelect_Parallel(y, dm, method = c("parallel"), mem.efficient = FALSE,
#'                   num_workers, alpha, corMat = NA, corMethod = "standard", max_ord=5, verbose = FALSE,
#'                   directed = FALSE)
PCSelect_Parallel <- function(y, dm, method = c("parallel"), mem.efficient = FALSE,
                               num_workers, alpha, corMat = NA, corMethod = "standard",max_ord=5, verbose = FALSE,
                               directed = FALSE){
   # print("----------in the pc_select_parallel myself's----------")
  time0 <- Sys.time()
  stopifnot((n <- nrow(dm)) >= 1, (p <- ncol(dm)) >= 1)
  vNms <- colnames(dm)
  # print("vNms = ")
  # print(vNms)
  cl <- match.call()
  zMin <- c(0, rep.int(Inf, p))
  #********************************************************************
  # Revised by:yinly
  # Date: Feb 27, 2020
  # Description: use a faster function to calculate correlation matrix
  #********************************************************************
  #C <- pcalg::mcor(cbind(y, dm), method = corMethod)
  #if (!is.na(corMat)) {C <- corMat} else {
  if (length(is.na(corMat)==0)){
    C <- corMat
  } else {
    C <- coop::pcor(cbind(y,dm))
  }
  
  cutoff <- qnorm(1 - alpha/2)
  n.edgetests <- numeric(1)
  G <- c(FALSE, rep.int(TRUE, p))
  seq_p <- 1:(p + 1)
  done <- FALSE
  ord <- 0
  if (num_workers < 2) {
    stop("The number of cores is insufficient to run parallel-PC")
  }
  workers <- NULL
  time1 <- Sys.time()
  if (Sys.info()[["sysname"]] == "Windows"){
    workers <- makeCluster(num_workers, type = "PSOCK")
    eval(dm)
    clusterEvalQ(workers, library(pcalg))
  }
  #print("parallel prepare cost time:")
  #print(Sys.time() - time1)
  time2 <- Sys.time()
  #print("parallel perpare finished.")
  edge_test_xy <- function(x, y){
    if (verbose >= 1) {
      #cat("in_fun edge_test_xy,x=", x, "y=", y)
    }
    G_x <- TRUE
    G_y <- TRUE
    num_tests_xy <- 0
    zMin_x <- zMin[x]
    done_xy <- TRUE
    if (G_x){
      nbrsBool <- as.logical(G.l)
      nbrsBool[x] <- FALSE
      nbrs <- seq_p[nbrsBool]
      length_nbrs <- length(nbrs)
      if (length_nbrs >= ord) {
        if (length_nbrs > ord)
          done_xy <- FALSE
        S <- seq(length = ord)
        repeat {
          num_tests_xy <- num_tests_xy + 1
          z <- pcalg::zStat(x, y, nbrs[S], C, n)
          if (abs(z) < zMin_x)
            zMin_x <- abs(z)
          if (verbose >= 2) {
            #print(paste0("here,x = ", x))
            #print(paste0("here,y = ", y))
            #cat(paste("x:", vNms[x - 1], "y:", (ytmp <- round((p + 1)/2)), "S:"), c(ytmp, vNms)[nbrs[S]],paste("z:", z, "\n"))
          }
          if (abs(z) <= cutoff) {
            G_x <- FALSE
            break
          }
          else {
            nextSet <- getNextSet(length_nbrs, ord, S)
            if (nextSet$wasLast) {
              break
            }
            S <- nextSet$nextSet
          }
        }
      }
    }
    list(G_x, num_tests_xy, zMin_x, done_xy)
  }
  edge_test <- function(i) {
    y <- 1
    x <- ind[i]
    if (verbose >= 1) {
      #cat("-----in edge test no xy-----i, x =", i, x, "\n")
    }
    num_tests_i <- 0
    G_i <- TRUE
    zMin_x <- zMin[x]
    done_i <- TRUE
    if (verbose >= 1) {
      #print("-----in edge_test(x, y)-----")
    }
    res_x <- edge_test_xy(x, y)
    G_i <- res_x[[1]]
    num_tests_i <- num_tests_i + res_x[[2]]
    zMin_x <- res_x[[3]]
    done_i <- done_i & res_x[[4]]
    list(i, G_i, num_tests_i, zMin_x, done_i)
  }
  edge_tests <- function(l) {
    res <- vector("list", length(l))
    for (k in 1:length(l)) {
      res[[k]] <- edge_test(l[[k]])
    }
    res
  }
  total_mem <- function() {
    tryCatch({
      if (Sys.info()[["sysname"]] == "Linux") {
        total <- (as.numeric(gsub("[^0-9]*([0-9]*)",
                                  "\\1", system("egrep '^MemFree:' /proc/meminfo",
                                                intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)",
                                                                                   "\\1", system("egrep '^Cached:' /proc/meminfo",
                                                                                                 intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)",
                                                                                                                                    "\\1", system("egrep '^Inactive:' /proc/meminfo",
                                                                                                                                                  intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)",
                                                                                                                                                                                     "\\1", system("egrep '^Buffers:' /proc/meminfo",
                                                                                                                                                                                                   intern = TRUE))))/1000
        return(total)
      }
      else if (Sys.info()[["sysname"]] == "Windows") {
        total <- (as.numeric(gsub("[^0-9]*([0-9]*)",
                                  "\\1", system("wmic OS get FreePhysicalMemory /Value",
                                                intern = TRUE))[3]))/1000
        return(total)
      }
      else if (Sys.info()[["sysname"]] == "Darwin") {
        total <- 4096 * (as.numeric(gsub("[^0-9]*([0-9]*)",
                                         "\\1", system("vm_stat | grep 'Pages free'",
                                                       intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)",
                                                                                          "\\1", system("vm_stat | grep 'Pages inactive'",
                                                                                                        intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)",
                                                                                                                                           "\\1", system("vm_stat | grep 'Pages speculative'",
                                                                                                                                                         intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)",
                                                                                                                                                                                            "\\1", system("vm_stat | grep 'Pages purgeable'",
                                                                                                                                                                                                          intern = TRUE))))/1e+06
        return(total)
      }
      else {
        total <- (as.numeric(gsub("[^0-9]*([0-9]*)",
                                  "\\1", system("vmstat -s | grep 'free memory'",
                                                intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)",
                                                                                   "\\1", system("vmstat -s | grep 'inactive memory'",
                                                                                                 intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)",
                                                                                                                                    "\\1", system("vmstat -s | grep 'buffer memory'",
                                                                                                                                                  intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)",
                                                                                                                                                                                     "\\1", system("vmstat -s | grep 'swap cache'",
                                                                                                                                                                                                   intern = TRUE))))/1000
        return(total)
      }
    }, error = function(e) {
      return(1024)
    }, warning = function(e) {
      return(1024)
    })
  }
  #print("-----pepare to parallel------")
  parallel_threshold <- 4
  if (mem.efficient) {
    mem_per_test <- 2
    tests_per_batch <- as.integer(total_mem()/mem_per_test)
  }
  time3 <- Sys.time()
  start_total <- proc.time()
  # ******************************************************
  # Revised by:yinly, Date: Feb 27, 2020
  # Description: use max_ord to control the workflow
  # ******************************************************
  while (!done && any(G) && ord<=max_ord) {
    ##G_order_specific stores which nodes are connected to the outcome for each order K
    G_order_specific <- G    
    
    n.edgetests[ord + 1] <- 0
    done <- TRUE
    # **************************************************
    # Revised by;yinly, Date: Feb 27, 2020
    # Description: replace G with G_order_specific
    # **************************************************
    ind <- which(G_order_specific) # which(G)
    remainingEdgeTests <- length(ind)
    if (verbose >= 1)
      cat("Order=", ord, "; remaining edges:", remainingEdgeTests,
          "\n", sep = "")
    # **************************************************
    # Revised by;yinly, Date: Feb 27, 2020
    # Description: replace G with G_order_specific
    # **************************************************
    G.l <- split(G_order_specific, gl(p + 1, 1)) # split(G, gl(p + 1, 1)
    if (!mem.efficient) {
      tests_per_batch <- remainingEdgeTests
    }
    for (j in seq(1, remainingEdgeTests, by = tests_per_batch)) {
      l <- min(remainingEdgeTests, j + tests_per_batch -
                 1)
      if (l - j + 1 < num_workers) {
        num_workers <- l - j + 1
      }
      res <- NULL
      if (l - j + 1 < parallel_threshold) {
        #cat("l-j+1 = :", l - j + 1, "parallel_threshold is :",parallel_threshold)
        #print("*** system ***")
        res <- lapply(j:l, edge_test)
      }
      else if (Sys.info()[["sysname"]] == "Windows") {
        #print("windows")
        res <- do.call("c", clusterApply(workers, clusterSplit(workers,
                                                               j:l), edge_tests))
      }
      else {
        #print("not windows")
        res <- mclapply(j:l, edge_test, mc.cores = num_workers,
                        mc.set.seed = FALSE, mc.cleanup = TRUE, mc.allow.recursive = FALSE)
      }
      for (p_obj in res) {
        i <- p_obj[[1]]
        y <- 1
        x <- ind[i]
        n.edgetests[ord + 1] <- n.edgetests[ord + 1] +
          p_obj[[3]]
        zMin[x] <- p_obj[[4]]
        G[x] <- p_obj[[2]]
        done <- done & p_obj[[5]]
      }
    }
    ord <- ord + 1
  }
  #print("parallel calculate cost:")
  #print(Sys.time() - time3)
  total_t = proc.time() - start_total
  #cat("Num CI Tests=", n.edgetests, ",Total CI Tests=", sum(unlist(n.edgetests)),",Total Time=", total_t[3], "\n", sep = " ")
  #cat("total is :", total_t, "\n")
  #cat("cost time 1:", time1 - time0, "cost time 2:", time2 - time1, "\n")
  Gres <- G[-1]
  names(Gres) <- vNms
  list(G = Gres, zMin = zMin[-1])
}
