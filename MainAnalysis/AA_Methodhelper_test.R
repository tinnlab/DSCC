### use euc_dist, and ang_dist to generate two combined affinity matrices, run spectral clustering on each to get one cluster
### also use aasc_eigengap to combine the affinity matrices, and run spectral clustering on the combined one
### add some parameters to decide whether using all metabolites or only ones can be mapped to pathway
### select best cluster among the three using silhouette score

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, NUMEXPR_NUM_THREADS = 1)

# Wrapper to run each method
suppressMessages({ suppressWarnings({
  library(wordspace)
  library(tidyverse)
  library(survival)
  library(matrixStats)
  library(SNFtool)
  library(igraph)
  library(NEMO)
  library(psych)
  library(parallel)
  library(cluster)
  library(philentropy)
}) })

clust_to_conn <- function(clustering) {
  n <- length(clustering)
  connectivity_matrix <- matrix(0, nrow = n, ncol = n)

  unique_clusters <- unique(clustering)
  for (cluster in unique_clusters) {
    members <- which(clustering == cluster)
    connectivity_matrix[members, members] <- 1
  }
  return(connectivity_matrix)
}

angular_distance_matrix <- function(data) {
  norms <- sqrt(rowSums(data^2))
  if (any(norms == 0)) {
    stop("Cannot compute angular distance with zero vector")
  }
  normalized_data <- data / norms
  cos_sim_matrix <- normalized_data %*% t(normalized_data)
  cos_sim_matrix <- pmin(pmax(cos_sim_matrix, -1.0), 1.0)
  angle_matrix <- acos(cos_sim_matrix)
  dist_matrix <- angle_matrix / pi
  diag(dist_matrix) <- 0
  return(dist_matrix)
}

correlation_distance_matrix <- function(data) {
  matrix_data <- t(data)
  cor_matrix <- cor(matrix_data)
  dist_matrix <- 1 - cor_matrix
  return(dist_matrix)
}

distance_to_affinity <- function(dist_matrix, sigma = 1) {
  affinity_matrix <- exp(-(dist_matrix^2) / (2 * sigma^2))
  return(affinity_matrix)
}

merge_matrices <- function(matrices) {
  n <- nrow(matrices[[1]])
  k <- length(matrices)
  w <- array(0, dim = c(n, n, k))
  for (i in 1:k) w[, , i] <- matrices[[i]]
  return(w)
}

solve_lambda <- function(alpha, gamma, ini) {

  f <- function(x) {
    a <- rep(1, length(gamma))
    abs((sum(a / ((gamma - x) * alpha))^2) /
          (sum(a / ((gamma - x)^2 * alpha))) - 1)
  }

  k <- length(gamma)
  # Initial x0
  x0 <- if (ini == 0) -1
  else if (ini == k + 1) 1
  else round(gamma[1] * 100) / 100
  # Use optim instead of fsolve (R equivalent)
  set.seed(123)
  opt_result <- optim(x0, f, method = "Brent",
                      lower = -2, upper = 2)
  return(opt_result$par)
}

aasc_eigengap <- function(w, num_clusters = 10) {
  sizeK <- dim(w)[3]  # Number of affinity matrices
  s1 <- dim(w)[1]     # Number of nodes
  s2 <- dim(w)[2]     # Should be equal to s1

  num_clusters <- min(num_clusters, s1)

  # Initialize matrices
  D_k <- array(0, dim = c(s1, s2, sizeK))
  L_k <- array(0, dim = c(s1, s2, sizeK))

  # Initial equal weight
  weight <- rep(1 / sizeK, sizeK)

  # Precompute degree and Laplacian matrices for each input matrix
  for (i in 1:sizeK) {
    # Compute degree matrix
    s <- rowSums(w[, , i])
    D_k[, , i] <- diag(s)
    # Compute Laplacian matrix
    L_k[, , i] <- D_k[, , i] - w[, , i]
  }

  # Parameters
  iter_number <- 20
  best_eigengap <- -Inf  # Changed from best_eigensum to best_eigengap
  best_weights <- weight

  # Main iterative process
  for (iter in 1:iter_number) {
    # Update weighted adjacency matrix
    # w_n <- matrix(0, nrow = s1, ncol = s2)
    # for (i in 1:sizeK) {
    #   w_n <- w_n + w[, , i] * (weight[i]^2)
    # }
    # w_n <- apply(w * (weight^2), MARGIN = c(1, 2), FUN = sum)

    w_reshaped <- matrix(w, nrow = s1 * s2, ncol = sizeK)
    # Multiply by weights and sum
    # w_n_flat <- w_reshaped %*% (weight^2)
    w_n_flat <- w_reshaped %*% weight
    # Reshape back to original dimensions
    w_n <- matrix(w_n_flat, nrow = s1, ncol = s2)

    # Compute degree matrix
    D <- diag(rowSums(w_n))

    # Compute normalized Laplacian matrix
    D_inv_sqrt <- diag(1 / sqrt(diag(D) + 1e-10))  # Add small constant to avoid division by zero
    L <- D - w_n
    L_norm <- D_inv_sqrt %*% L %*% D_inv_sqrt

    # Compute eigendecomposition
    eigen_result <- eigen(L_norm, symmetric = TRUE)
    eigenvalues <- eigen_result$values
    eigenvectors <- eigen_result$vectors

    # Sort eigenvalues and eigenvectors
    idx <- order(eigenvalues)
    eigenvalues <- eigenvalues[idx]
    eigenvectors <- eigenvectors[, idx]

    # Calculate all eigengaps (differences between consecutive eigenvalues)
    # We'll look at all possible gaps to find the most significant one
    eigengaps <- abs(diff(eigenvalues))
    eigengaps <- (1:length(eigengaps)) * eigengaps

    # For the purposes of clustering, we're most interested in the largest gap
    # among the first several eigenvalues (up to num_clusters*2 to be safe)
    # search_range <- min(length(eigengaps), 2*num_clusters)
    # max_gap_idx <- which.max(eigengaps[1:length(eigengaps)])
    eigengap <- max(eigengaps)

    # Keep track of the best combination found - now maximizing eigengap
    if (eigengap > best_eigengap) {
      best_eigengap <- eigengap
      best_weights <- weight
    }

    # Compute coefficients for updating weights
    # alpha_k <- matrix(0, nrow = num_clusters, ncol = sizeK)
    # beta_k <- matrix(0, nrow = num_clusters, ncol = sizeK)
    # gamma_k <- matrix(0, nrow = num_clusters, ncol = sizeK)

    relevant_eigenvectors <- eigenvectors[, 2:num_clusters, drop = FALSE]

    # for (k in 1:sizeK) {
    #   # For each eigenvector and each input matrix
    #   for (c in 1:num_clusters) {
    #     eigen_idx <- c + 1  # Skip the trivial eigenvector
    #     eigen_vec <- eigenvectors[, eigen_idx]
    #
    #     # Compute Rayleigh quotient
    #     alpha_k[c, k] <- t(eigen_vec) %*% D_k[, , k] %*% eigen_vec
    #     beta_k[c, k] <- t(eigen_vec) %*% L_k[, , k] %*% eigen_vec
    #     gamma_k[c, k] <- beta_k[c, k] / alpha_k[c, k]
    #   }
    # }

    alpha_k <- sapply(1:sizeK, function(k) {
      D_mult <- D_k[, , k] %*% relevant_eigenvectors
      colSums(relevant_eigenvectors * D_mult)
    })

    # Calculate beta_k values using sapply
    beta_k <- sapply(1:sizeK, function(k) {
      L_mult <- L_k[, , k] %*% relevant_eigenvectors
      colSums(relevant_eigenvectors * L_mult)
    })

    # Calculate gamma_k values in one operation
    gamma_k <- beta_k / alpha_k

    # Average across all eigenvectors of interest to get a combined measure
    alpha <- colMeans(alpha_k)
    gamma <- colMeans(gamma_k)

    # 1-D solution for weight optimization
    weight_1D <- matrix(0, nrow = sizeK + 2, ncol = sizeK)
    result_1D <- numeric(sizeK + 2)

    for (iter_1 in 0:(sizeK + 1)) {
      lambda_1 <- solve_lambda(alpha, gamma, iter_1)
      a <- rep(1, sizeK)
      b <- (gamma - lambda_1) * alpha
      lambda_2 <- 1 / sum(a / b)
      c <- sqrt(alpha) * (gamma - lambda_1)
      u <- a * lambda_2 / c
      weight_1D[iter_1 + 1,] <- u / sqrt(alpha)

      # Evaluate the quality of this weight combination - now using eigengap
      # test_w_n <- matrix(0, nrow = s1, ncol = s2)
      # for (i in 1:sizeK) {
      #   test_w_n <- test_w_n + w[, , i] * (weight_1D[iter_1 + 1, i]^2)
      # }

      w_reshaped <- matrix(w, nrow = s1 * s2, ncol = sizeK)
      test_w_n_flat <- w_reshaped %*% (weight_1D[iter_1 + 1,]^2)
      test_w_n <- matrix(test_w_n_flat, nrow = s1, ncol = s2)

      test_D <- diag(rowSums(test_w_n))
      test_L <- test_D - test_w_n
      test_D_inv_sqrt <- diag(1 / sqrt(diag(test_D) + 1e-10))
      test_L_norm <- test_D_inv_sqrt %*% test_L %*% test_D_inv_sqrt
      test_eigen <- eigen(test_L_norm, symmetric = TRUE)
      test_eigenvalues <- sort(test_eigen$values)

      # Calculate all eigengaps and find the maximum
      test_eigengaps <- abs(diff(test_eigenvalues))
      test_eigengaps <- (1:length(test_eigengaps)) * test_eigengaps
      # search_range <- min(length(test_eigengaps), 2*num_clusters)
      # result_1D[iter_1 + 1] <- max(test_eigengaps[1:search_range])
      result_1D[iter_1 + 1] <- max(test_eigengaps)

    }

    # Select weights that give maximum eigengap (not minimum eigensum)
    max_idx <- which.max(result_1D)
    weight <- weight_1D[max_idx,]

    # Normalize weights to sum to 1
    weight <- weight / sum(weight)
  }

  # Return the weights that gave the best result
  return(list(
    weights = best_weights,
    eigengap = best_eigengap # Return eigengap instead of eigensum
  ))
}


fuse_conn_matrices <- function(allSamples, conn.per.omic) {
  patient.names <- as.vector(allSamples)
  num.patients <- length(patient.names)
  returned.conn.matrix <- matrix(0, ncol = num.patients, nrow = num.patients)
  rownames(returned.conn.matrix) <- patient.names
  colnames(returned.conn.matrix) <- patient.names

  shared.omic.count <- matrix(0, ncol = num.patients, nrow = num.patients)
  rownames(shared.omic.count) <- patient.names
  colnames(shared.omic.count) <- patient.names

  for (j in 1:length(conn.per.omic)) {
    curr.omic.patients <- colnames(conn.per.omic[[j]])
    returned.conn.matrix[curr.omic.patients, curr.omic.patients] <- returned.conn.matrix[curr.omic.patients, curr.omic.patients] + conn.per.omic[[j]][curr.omic.patients, curr.omic.patients]
    shared.omic.count[curr.omic.patients, curr.omic.patients] <- shared.omic.count[curr.omic.patients, curr.omic.patients] + 1
  }

  final.ret <- returned.conn.matrix / shared.omic.count
  lower.tri.ret <- final.ret[lower.tri(final.ret)]
  final.ret[shared.omic.count == 0] <- mean(lower.tri.ret[!is.na(lower.tri.ret)])
  ## prevent the connectivity matrix from having NA values
  final.ret[is.na(final.ret)] <- 0
  diag(final.ret) <- 1
  final.ret
}

filterGenesUsingGeneList <- function(M, geneList, delimeter) {

  if (is.null(delimeter)) {
    features <- colnames(M)
  } else {
    # features <- colnames(M) %>% stringr::str_remove(paste0(delimeter, ".+$")) %>% toupper()
    features <- colnames(M) %>%
      lapply(function(x) {
        stringr::str_split(x, delimeter)[[1]][1]
      }) %>%
      unlist() %>%
      toupper()
  }

  if (length(features) %in% c(0, 1)) {
    return(NULL)
  }

  keepFeatures <- features %in% (geneList %>% unique())

  M <- as.matrix(M[, keepFeatures]) %>%
    as.data.frame() %>%
    t() %>%
    unique() %>%
    t()

  rn <- rownames(M)
  cn <- colnames(M)

  M <- M %>% as.numeric() %>% matrix(nrow = nrow(M))
  colnames(M) <- cn
  rownames(M) <- rn
  M
}

runDSCC <- function(dataList, nSamples, defk = 10) {

  dataList <- lapply(names(dataList), function(dataType) {
    data <- dataList[[dataType]]

    data <- as.data.frame(data) %>% t() %>% unique() %>% t()
    filteredDat <- data[, colVars(data) > 0, drop = FALSE]
    print(ncol(filteredDat))
    # filteredDat <- filteredDat[, order(colVars(filteredDat), decreasing = TRUE)[1:min(2000, ncol(filteredDat))]]
    # scale(filteredDat)
    # if(dataType != "meth450"){
    #   return(scale(filteredDat))
    # }else{
    #   return(filteredDat)
    # }
    filteredDat
  })

  ang_aff <- lapply(dataList, function(data) {
    dist_matrix <- angular_distance_matrix(as.matrix(data))
    # non_zero_dists <- dist_matrix[dist_matrix > 0]
    # aff <- distance_to_affinity(dist_matrix, 0.2)
    aff <- distance_to_affinity(dist_matrix, 0.5)
    # aff <- distance_to_affinity(dist_matrix, median(non_zero_dists))

    if (nrow(dist_matrix) <= defk){
      tmpk <- nrow(dist_matrix)
    }else{
      tmpk <- defk
    }
    if (nSamples > 200) {
      # non_zero_dists <- dist_matrix[dist_matrix > 0]
      # if (length(non_zero_dists) > 0) {
      #   sigma <- mean(non_zero_dists) / 1  # A reasonable heuristic
      # } else {
      #   sigma <- 1.0  # Fallback
      # }
    # # }else {
    #   # aff <- NEMO::affinityMatrix(dist_matrix, defk, 0.5)
      non.sym.knn <- apply(aff, 1, function(sim.row) {
        returned.row <- sim.row
        threshold <- sort(sim.row, decreasing = T)[tmpk]
        returned.row[sim.row < threshold] <- 0
        row.sum <- sum(returned.row)
        returned.row[sim.row >= threshold] <- returned.row[sim.row >= threshold] / row.sum
        # returned.row <- returned.row / row.sum
        return(returned.row)
      })
      aff <- non.sym.knn + t(non.sym.knn)
    }
    aff
  })

  mag_aff <- lapply(dataList, function(data) {
    dist_matrix <- as.matrix(dist(as.matrix(data), method = "euclidean"))
    if (nrow(dist_matrix) <= defk){
      tmpk <- nrow(dist_matrix)
    }else{
      tmpk <- defk
    }
    non_zero_dists <- dist_matrix[dist_matrix > 0]
    if (length(non_zero_dists) > 0) {
      sigma <- median(non_zero_dists) / 1  # A reasonable heuristic
    } else {
      sigma <- 1.0  # Fallback
    }

    # aff <- SNFtool::affinityMatrix(dist_matrix, tmpk, 0.5)

    aff <- distance_to_affinity(dist_matrix, sigma)
    if (nSamples > 200) {
      # aff <- distance_to_affinity(dist_matrix, sigma)
    # }else {
    #   aff <- SNFtool::affinityMatrix(dist_matrix, defk, sigma)
      non.sym.knn <- apply(aff, 1, function(sim.row) {
        returned.row <- sim.row
        threshold <- sort(sim.row, decreasing = T)[tmpk]
        returned.row[sim.row < threshold] <- 0
        row.sum <- sum(returned.row)
        returned.row[sim.row >= threshold] <- returned.row[sim.row >= threshold] / row.sum
        # returned.row <- returned.row / row.sum
        return(returned.row)
      })
      aff <- non.sym.knn + t(non.sym.knn)
    }
    aff
  })

  # cor_aff <- lapply(dataList, function(data) {
  #   dist_matrix <- correlation_distance_matrix(data)
  #   non_zero_dists <- dist_matrix[dist_matrix > 0]
  #   if (length(non_zero_dists) > 0) {
  #     sigma <- median(non_zero_dists) / 1  # A reasonable heuristic
  #   } else {
  #     sigma <- 1.0  # Fallback
  #   }
  #   distance_to_affinity(dist_matrix, sigma)
  #   # distance_to_affinity(dist_matrix, 0.5)
  # })

  aasc_res <- lapply(1:length(mag_aff), function(i) {
    # array3D <- merge_matrices(list(mag_aff[[i]], cor_aff[[i]], ang_aff[[i]]))
    array3D <- merge_matrices(list(mag_aff[[i]], ang_aff[[i]]))
    opt_weights <- aasc_eigengap(array3D, num_clusters = defk)$weights
    # aff <- mag_aff[[i]] * opt_weights[1] + cor_aff[[i]] * opt_weights[2] + ang_aff[[i]] * opt_weights[3]
    aff <- mag_aff[[i]] * (opt_weights[1]) + ang_aff[[i]] * (opt_weights[2])
    rownames(aff) <- colnames(aff) <- rownames(mag_aff[[i]])

    num.clusters <- nemo.num.clusters(aff)
    # num.clusters <- 4
    cluster <- spectralClustering(aff, num.clusters)
    names(cluster) <- rownames(aff)
    conn <- clust_to_conn(cluster)
    rownames(conn) <- colnames(conn) <- rownames(aff)
    return(list(aff = aff, conn = conn))
  })

  allSamples <- lapply(mag_aff, function(sim) rownames(sim)) %>%
    unlist() %>%
    unique()

  affinity.per.omic <- lapply(aasc_res, function(res) { res$aff })
  conn.per.omic <- lapply(aasc_res, function(res) { res$conn })

  final.mag.aff <- fuse_conn_matrices(allSamples, mag_aff)
  # final.cor.aff <- fuse_conn_matrices(allSamples, cor_aff)
  final.ang.aff <- fuse_conn_matrices(allSamples, ang_aff)
  # final.aff.dist <- (final.mag.aff + final.ang.aff)/2
  final.aff <- fuse_conn_matrices(allSamples, affinity.per.omic)
  final.conn <- fuse_conn_matrices(allSamples, conn.per.omic)

  ## Cluster 1
  num.clusters <- nemo.num.clusters(final.mag.aff, 2:10)
  # num.clusters <- 4
  cluster1 <- spectralClustering(final.mag.aff, num.clusters)
  names(cluster1) <- allSamples
  dist_mat1 <- 1 - final.aff
  sh_1 <- try(mean(silhouette(cluster1, dist_mat1)[, "sil_width"]), silent = TRUE)

  ## Cluster 2
  # num.clusters <- nemo.num.clusters(final.cor.aff)
  # # num.clusters <- 4
  # cluster2 <- spectralClustering(final.cor.aff, num.clusters)
  # names(cluster2) <- allSamples
  # dist_mat2 <- 1-final.aff
  # sh_2 <- try(mean(silhouette(cluster2, dist_mat2)[, "sil_width"]), silent = TRUE)

  ## Cluster 3
  num.clusters <- nemo.num.clusters(final.ang.aff, 2:10)
  # num.clusters <- 4
  cluster3 <- spectralClustering(final.ang.aff, num.clusters)
  names(cluster3) <- allSamples
  dist_mat3 <- 1 - final.aff
  sh_3 <- try(mean(silhouette(cluster3, dist_mat3)[, "sil_width"]), silent = TRUE)

  ## Cluster 4
  # num.clusters <- nemo.num.clusters(final.aff)
  # # num.clusters <- 4
  # cluster4 <- spectralClustering(final.aff, num.clusters)
  # names(cluster4) <- allSamples
  # dist_mat4 <- 1-final.aff.dist
  # sh_4 <- try(mean(silhouette(cluster4, dist_mat4)[, "sil_width"]), silent = TRUE)

  ## Cluster 5
  probs <- seq(0.01, 0.99, by = 0.01) # Creates sequence from 1% to 99%
  quantiles <- quantile(final.conn, probs = probs)
  mod_optimum <- 0
  dist_mat5 <- 1 - final.aff

  for (connThreshold in quantiles) {
    conn <- final.conn
    conn[conn < connThreshold] <- 0
    # dist_mat_tmp <- as.dist(1-conn)
    graph <- graph_from_adjacency_matrix(conn, mode = "undirected", weighted = T)
    for (seed in 1:50) {
      set.seed(seed)
      louvain_res <- cluster_louvain(graph)
      cluster <- membership(louvain_res)  # Use membership() function
      names(cluster) <- allSamples
      # print(max(cluster))
      # Calculate silhouette score
      # current_sh <- try(mean(silhouette(cluster, dist_mat)[, "sil_width"]), silent = TRUE)
      current_modularity <- modularity(graph, cluster)
      if (!inherits(current_modularity, "try-error")) {
        if (max(cluster) <= 10 & current_modularity > mod_optimum) {
          # if (current_sh > sh_optimum){
          mod_optimum <- current_modularity
          cluster5 <- cluster
          # dist_mat2 <- dist_mat_tmp
        }
      }
    }
  }
  sh_5 <- try(mean(silhouette(cluster5, dist_mat5)[, "sil_width"]), silent = TRUE)
  if (inherits(sh_5, "try-error")) {
    sh_5 <- 0
    cluster5 <- NULL
  }

  # ## Cluster 3
  # num.clusters <- nemo.num.clusters(final.cor.aff)
  # # num.clusters <- 4
  # cluster3 <- spectralClustering(final.cor.aff, num.clusters)
  # names(cluster3) <- allSamples
  # dist_mat3 <- 1-final.aff
  # sh_3 <- try(mean(silhouette(cluster3, dist_mat3)[, "sil_width"]), silent = TRUE)
  #
  all_sh <- c(sh_1, sh_3, sh_5)
  all_cluster <- list(cluster1, cluster3, cluster5)

  return(all_cluster[[ which.max(all_sh)[1] ]])

  # return(cluster1)
  # if (sh_1 >= sh_2){
  #   return(cluster1)
  # }else{
  #   return(cluster1)
  # }
}

