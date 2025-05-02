# Wrapper to run each method
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, NUMEXPR_NUM_THREADS = 1)

# Wrapper to run each method
suppressMessages({ suppressWarnings({
  library(ANF)
  library(PMA)
  library(IntNMF)
  library(wordspace)
  library(kernlab)
  library(tidyverse)
  library(survival)
  library(matrixStats)
  library(ConsensusClusterPlus)
  library(CIMLR)
  library(SNFtool)
  library(iClusterPlus)
  library(NNLM)
  library(igraph)
  library(NEMO)
  library(PINSPlus)
  library(mogsa)
  library(polycor)
  library(psych)
  library(FactoMineR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pbmcapply)
  library(LRAcluster)
  library(cluster)
  library(survminer)
  # library(MOVICS)

}) })

LRAClus_prob <- readRDS("/data/daotran/Cancer_RP/Subtyping/Data/LRAClus_meth.rds")

# BiocManager::install(c("iClusterPlus", "mogsa"))
# mamba install -c conda-forge r-igraph r-polycor r-pbmcapply r-base=4.1.3
# install.packages(c("PINSPlus", "psych", "FactoMineR"))

standardNormalization <- function(x) {
  x = as.matrix(x)
  mean = apply(x, 2, mean)
  sd = apply(x, 2, sd)
  sd[sd == 0] = 1
  xNorm = t((t(x) - mean) / sd)
  return(xNorm)
}

listNorm <- function(dataList){
  newdataList <- lapply(names(dataList), function(dataType){
    data <- dataList[[dataType]]
    if(!dataType %in% c("meth450", "clinicalMCA")){
      # if(dataType == "clinicalImputedV2"){
      #   # num_cols <- c("age_at_initial_pathologic_diagnosis", "height", "weight")
      #   # keeps <- lapply(colnames(data), function(name){
      #   #   countmatch <- sapply(num_cols, function(col) {
      #   #     length(grep(col, name))
      #   #   }) %>% sum()
      #   #   if(countmatch == 0){
      #   #     return(NULL)
      #   #   }else{
      #   #     return(name)
      #   #   }
      #   # }) %>% do.call(what = c) %>% unique()
      #   # keeps <- keeps[!is.null(keeps)]
      #   # if (length(keeps) > 0){
      #   #   for (i in keeps){
      #   #     data[, i] <- scale(data[, i])
      #   #   }
      #   # }
      #   for (i in 1:ncol(data)){
      #     if (max(data[, i]) > 1){
      #       data[, i] <- scale(data[, i])
      #     }
      #   }
      # }else{
        if(min(data) >= 0){
          data <- standardNormalization(data)
        }
      # }
    }
    return(as.matrix(data))
  }) %>% `names<-` (names(dataList))
  return(newdataList)
}


#### CC ####
runCC <- function(dataList) {

  triangle = function(m, mode = 1) {
    #mode=1 for CDF, vector of lower triangle.
    #mode==3 for full matrix.
    #mode==2 for calcICL; nonredundant half matrix coun
    #mode!=1 for summary
    n = dim(m)[1]
    nm = matrix(0, ncol = n, nrow = n)
    fm = m


    nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half

    fm = t(nm) + nm
    diag(fm) = diag(m)

    nm = fm
    nm[upper.tri(nm)] = NA
    diag(nm) = NA
    vm = m[lower.tri(nm)]

    if (mode == 1) {
      return(vm) #vector
    }else if (mode == 3) {
      return(fm) #return full matrix
    }else if (mode == 2) {
      return(nm) #returns lower triangle and no diagonal. no double counts.
    }
  }

  maxK <- 10

  # d=cbind(dataList[[1]], dataList[[2]], dataList[[3]])


  dataList <- listNorm(dataList)
  d = do.call(cbind, dataList)

  d = standardNormalization(t(d))
  # d = standardNormalization(d); d = t(d)
  # d = sweep(d, 1, apply(d, 1, median, na.rm = T))
  set.seed(1234)
  result_All = ConsensusClusterPlus(d, maxK = maxK, reps = 50, pItem = 0.8, pFeature = 1, clusterAlg = "hc", distance = "pearson", plot = NULL, seed = 1)

  k = maxK
  breaks = 100
  areaK = c()
  for (i in 2:k) {
    v = triangle(result_All[[i]]$ml, mode = 1)

    #empirical CDF distribution. default number of breaks is 100
    h = hist(v, plot = FALSE, breaks = seq(0, 1, by = 1 / breaks))
    h$counts = cumsum(h$counts) / sum(h$counts)

    #calculate area under CDF curve, by histogram method.
    thisArea = 0
    for (bi in 1:(length(h$breaks) - 1)) {
      thisArea = thisArea + h$counts[bi] * (h$breaks[bi + 1] - h$breaks[bi]) #increment by height by width
      bi = bi + 1
    }
    areaK = c(areaK, thisArea)
  }

  #plot area under CDF change.
  deltaK = areaK[1] #initial auc at k=2
  for (i in 2:(length(areaK))) {
    #proportional increase relative to prior K.
    deltaK = c(deltaK, (areaK[i] - areaK[i - 1]) / areaK[i - 1])
  }

  k <- which.max(deltaK[2:(k - 1)]) + 2 #minimun no cluster is 3, first pos is not calculated (k=2)
  cluster <- result_All[[k]]$consensusClass

  cluster
}

#### SNF ####
runSNF <- function(dataList) {
  ## First, set all the parameters:
  K = 10; ##number of neighbors, usually (10~30)
  alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
  NIT = 20; ###Number of Iterations, usually (10~20)

  # Normalization
  dataList <- listNorm(dataList)
  PS <- lapply(dataList, function(dat) dist2(as.matrix(dat), as.matrix(dat)))

  # Add code to handle the problem of low number of samples
  minSamp <- min(lapply(dataList, nrow) %>% do.call(what = c))
  if (minSamp <= 20) {
    K <- minSamp - 1
  }

  W <- lapply(PS, function(ps) affinityMatrix(ps, K, alpha))
  # W = SNF(W, K, NIT)
  if (length(W) > 1) {
    W = SNF(W, K, NIT)
  } else {
    W = W[[1]]
  }

  #####Clustering
  #Groups with SNF
  C = estimateNumberOfClustersGivenGraph(W, 2:10)[[1]]

  set.seed(1234)
  groupSNF = spectralClustering(W, C); #the final subtypes information
  names(groupSNF) <- rownames(W)
  groupSNF
}

runCIMLR <- function(dataList) {
  NUMC = 2:10
  dataList <- listNorm(dataList)
  dataList <- lapply(dataList, function(data) {
    data[data > 10] <- 10
    data[data < -10] <- -10
    t(data)
  })

  set.seed(1234)
  k <- CIMLR_Estimate_Number_of_Clusters(dataList, NUMC, cores.ratio = 0)
  k <- NUMC[which.min(k$K1)]

  # closeAllConnections()

  set.seed(1234)
  result <- CIMLR(dataList, k, cores.ratio = 0)
  cluster <- result$y$cluster
  names(cluster) <- colnames(dataList[[1]])
  cluster
}

### NEMO ####
runNEMO <- function(dataList) {

  nemo.affinity.graph <- function(raw.data, k = NA) {
    if (is.na(k)) {
      k = as.numeric(lapply(1:length(raw.data), function(i) round(ncol(raw.data[[i]]) / NUM.NEIGHBORS.RATIO)))
    } else if (length(k) == 1) {
      k = rep(k, length(raw.data))
    }
    k[k <= 1] <- 2

    # k <- as.numeric(lapply(1:length(raw.data), function(i){
    #   return(min(ncol(raw.data[[i]])-1, 20))
    # }))
    # k[k <= 1] <- 2

    sim.data = lapply(1:length(raw.data), function(i) { affinityMatrix(dist2(as.matrix(t(raw.data[[i]])),
                                                                             as.matrix(t(raw.data[[i]]))), k[i], 0.5) })
    affinity.per.omic = lapply(1:length(raw.data), function(i) {
      sim.datum = sim.data[[i]]
      non.sym.knn = apply(sim.datum, 1, function(sim.row) {
        returned.row = sim.row
        threshold = sort(sim.row, decreasing = T)[k[i]]
        returned.row[sim.row < threshold] = 0
        row.sum = sum(returned.row)
        returned.row[sim.row >= threshold] = returned.row[sim.row >= threshold] / row.sum
        return(returned.row)
      })
      sym.knn = non.sym.knn + t(non.sym.knn)
      return(sym.knn)
    })
    patient.names = Reduce(union, lapply(raw.data, colnames))
    num.patients = length(patient.names)
    returned.affinity.matrix = matrix(0, ncol = num.patients, nrow = num.patients)
    rownames(returned.affinity.matrix) = patient.names
    colnames(returned.affinity.matrix) = patient.names

    shared.omic.count = matrix(0, ncol = num.patients, nrow = num.patients)
    rownames(shared.omic.count) = patient.names
    colnames(shared.omic.count) = patient.names

    for (j in 1:length(raw.data)) {
      curr.omic.patients = colnames(raw.data[[j]])
      returned.affinity.matrix[curr.omic.patients, curr.omic.patients] = returned.affinity.matrix[curr.omic.patients, curr.omic.patients] + affinity.per.omic[[j]][curr.omic.patients, curr.omic.patients]
      shared.omic.count[curr.omic.patients, curr.omic.patients] = shared.omic.count[curr.omic.patients, curr.omic.patients] + 1
    }

    final.ret = returned.affinity.matrix / shared.omic.count
    lower.tri.ret = final.ret[lower.tri(final.ret)]
    final.ret[shared.omic.count == 0] = mean(lower.tri.ret[!is.na(lower.tri.ret)])

    return(final.ret)
  }

  dataList <- listNorm(dataList)
  dataList <- lapply(dataList, function(data){
    t(data)
  })
  graph <- nemo.affinity.graph(dataList)
  num.clusters <- nemo.num.clusters(graph, 2:10)
  clustering <- spectralClustering(graph, num.clusters)
  names(clustering) <- colnames(graph)
  clustering
}

### LRACluster ###
runLRACluster <- function(dataList) {
  if ("meth450" %in% names(dataList)) {
    tmp <- dataList$meth450
    dataList$meth450 <- tmp[, colnames(tmp) %in% LRAClus_prob]
    rm(tmp)
  }

  dataList <- lapply(dataList, t)
  types = c()
  for (data in dataList) {
    if (all(data %in% 0:1)) {
      types = c(types, "binary")
    }
    else {
      types = c(types, "gaussian")
    }
  }

  rlist <- LRAcluster::LRAcluster(data = dataList, types = types)
  df <- t(rlist$coordinate)

  silhouette_score <- function(k) {
    set.seed(1234)
    km <- kmeans(df, centers = k)
    ss <- silhouette(km$cluster, dist(df))
    mean(ss[, 3])
  }

  k <- 2:10
  avg_sil <- sapply(k, silhouette_score)

  set.seed(1234)
  clst <- kmeans(df, which.max(avg_sil) + 1)
  clst$cluster
}

runIntNMF <- function(dataList) {
  # keeps <- lapply(names(dataList), function(name) {
  #   countmatch <- sapply(c("mRNA", "GeneExp"), function(dataType) {
  #     length(grep(dataType, name))
  #   }) %>% sum()
  #   if (countmatch == 0) {
  #     return(NULL)
  #   }else {
  #     return(name)
  #   }
  # }) %>%
  #   do.call(what = c) %>%
  #   unique()
  # keeps <- keeps[!is.null(keeps)]
  # if (length(keeps) > 0) {
  #   for (name in keeps) {
  #     tmp <- as.matrix(dataList[[name]])
  #     # tmp[is.na(tmp)] <- 0
  #     sds <- colSds(tmp)
  #     if (sum(sds >= 1.5) >= 2000){
  #       dataList[[name]] <- tmp[, sds >= 1.5]
  #     }else{
  #       dataList[[name]] <- tmp[, order(sds, decreasing = TRUE)[1:2000]]
  #     }
  #   }
  # }

  # if ("miRNA" %in% names(dataList)) {
  #   tmp <- as.matrix(dataList$miRNA)
  #   # tmp[is.na(tmp)] <- 0
  #   tmp <- tmp[, colMeans(tmp == 0) < 0.5]
  #   dataList$miRNA <- tmp
  #   rm(tmp)
  # }

  # if ("cnv" %in% names(dataList)) {
  #   tmp <- as.matrix(dataList$cnv)
  #   sds <- colSds(tmp)
  #   # tmp[is.na(tmp)] <- 0
  #   if (sum(sds >= 0.9) >= 2000){
  #     dataList$cnv <- tmp[, sds >= 0.9]
  #   }else{
  #     dataList$cnv <- tmp[, order(sds, decreasing = TRUE)[1:2000]]
  #   }
  #   rm(tmp, sds)
  # }

  # keep maximum 2000 genes
  maxGeneNo = 2000
  dataList <- lapply(dataList, function(dat) {
    dat <- as.matrix(dat)
    if (ncol(dat) > maxGeneNo) {
      sds <- colSds(dat)
      dat <- dat[, order(sds, decreasing = T)[1:maxGeneNo]]
    }
    dat
  })
  # making negative to positive
  dataList <- lapply(dataList, function(dat) {
    if (min(dat) < 0){
      dat <- (dat - min(dat)) # make it non negative
    }
    dat[, colSds(dat) > 0 & colSums(dat) > 0] + .Machine$double.eps
  })

  # Finding the optimum number of cluster
  set.seed(1234)
  kList <- IntNMF::nmf.opt.k(dataList, n.runs = 5, n.fold = 5, k.range = 2:10, make.plot = F)

  k <- which.max(rowMeans(kList)) + 1

  set.seed(1234)
  clusters <- IntNMF::nmf.mnnals(dataList, k = k)$clusters

  clusters
}

runANF <- function(dataList) {

  ## First, set all the parameters:
  K = 10; ##number of neighbors, usually (10~30)

  # Normalization
  dataList <- listNorm(dataList)
  PS <- lapply(dataList, function(dat) dist2(as.matrix(dat), as.matrix(dat)))
  WList <- lapply(PS, function(ps) ANF::affinity_matrix(ps, K))
  set.seed(1234)
  W <- ANF::ANF(WList, K = K)
  C = estimateNumberOfClustersGivenGraph(W, 2:10)[[1]]
  set.seed(1234)
  cluster = spectral_clustering(W, C); #the final subtypes information
  names(cluster) <- rownames(W)
  cluster

}