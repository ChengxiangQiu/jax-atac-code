### tf is peak * cell, idf is peak * 1, this function is used to let each column of tf multiply by the idf
### it means every column multiplys a same vector (idf)
safe_tfidf_multiply = function(tf, idf) {
  tf = t(tf)
  tf@x <- tf@x * rep.int(idf, diff(tf@p))
  tf = t(tf)
  return(tf)
}

## scale_vector should be a vector with length of column number, means each column of bmat divide by the the corresponding element of scale_vector
safe_column_scale = function(bmat, scale_vector) {
  bmat@x <- bmat@x / rep.int(scale_vector, diff(bmat@p)) 
  return(bmat)
}

safe_column_multiply = function(bmat, scale_vector) {
  bmat@x <- bmat@x * rep.int(scale_vector, diff(bmat@p))
  return(bmat)
}
# Perform TF-IDF on binary matrix
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
#   frequencies (bool): divide bmat by colSums (if FALSE simply use bmat for TF matrix)
#   log_scale_tf (bool): log scale TF matrix if TRUE
#   scale_factor (float): multiply terms in TF matrix by scale_factor prior to log1p. Equivalent to adding small pseudocount but doesn't cast to dense matrix at any point.
# Returns:
#   sparse matrix: TF-IDF matrix
tfidf = function(bmat, frequencies=TRUE, log_scale_tf=TRUE, scale_factor=100000) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  # TF-IDF
  tf_idf_counts = safe_tfidf_multiply(tf, idf)
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  return(tf_idf_counts)
}

sparse_prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, ...) {
  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  args <- list(A=x, nv=n, scale.=scale., center=center)
  if (!missing(...)) args <- c(args, list(...))
  s <- do.call(irlba, args=args)
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  ans$center <- args$center
  if (retx)
  {
    ans <- c(ans, list(d=s$d, x = s$u %*% diag(s$d)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}

atac_scrublet = function(bmat, k=NULL, fraction_sim_doublets=2, estimated_doublet_rate=0.1, dims=2:50) {
  # Determine KNN parameters
  if (is.null(k)) {
    k = round(0.5 * sqrt(ncol(bmat)))
  }
  kadj = round(k * (1+fraction_sim_doublets))
  # Perform TFIDF on original dataset
  message('[scrublet atac] Performing LSI-logTF on dataset...')
  # TF-IDF on original dataset
  bmat_colsums = Matrix::colSums(bmat)
  tf_original = safe_column_scale(bmat, bmat_colsums)
  tf_original@x = log1p(tf_original@x * 100000)
  idf_original = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  bmat.tfidf = safe_tfidf_multiply(tf_original, idf_original)
  rm(tf_original)
  bmat.pca = sparse_prcomp_irlba(t(bmat.tfidf), n=max(dims), center=FALSE, scale.=FALSE)
  rm(bmat.tfidf)
  # Make simulated doublets
  message('[scrublet atac] Simulating doublets...')
  set.seed(2019)
  doublet_half1 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
  doublet_half2 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
  doublets.bmat = bmat[, doublet_half1] + bmat[, doublet_half2]
  doublets.bmat@x[doublets.bmat@x > 1] = 1
  colnames(doublets.bmat) = paste0('DOUBLET_', 1:ncol(doublets.bmat))
  # Perform TF-IDF on doublets using IDF from original dataset
  doublet_colsums = bmat_colsums[doublet_half1] + bmat_colsums[doublet_half2] ## approximate (because of binarization after sum) to save time, but could recalculate
  tf_doublets = safe_column_scale(doublets.bmat, doublet_colsums)
  rm(doublets.bmat)
  tf_doublets@x = log1p(tf_doublets@x * 100000)
  doublets.tfidf = safe_tfidf_multiply(tf_doublets, idf_original)
  rm(tf_doublets)
  # Project doublets into PCA space and weight by variance explained
  message('[scrublet atac] Projecting doublets into PCA space...')
  doublets.pca = t(doublets.tfidf) %*% bmat.pca$rotation
  rm(doublets.tfidf)
  doublets.pca.weighted = doublets.pca
  # Jam it all into a Seurat object for some of the downstream stepst
  message('[scrublet atac] Making Seurat object...')
  rownames(bmat.pca$x) = colnames(bmat)
  combined.metadata = rbind(data.frame(cell=colnames(bmat), doublet=FALSE), data.frame(cell=rownames(doublets.pca), doublet=TRUE))
  rownames(combined.metadata) = combined.metadata$cell
  combined.pca = rbind(bmat.pca$x, doublets.pca.weighted)
  rownames(combined.pca) = rownames(combined.metadata)
  combined.seurat = CreateSeuratObject(counts=t(combined.pca), meta.data = combined.metadata)
  combined.seurat[['pca']] = Seurat::CreateDimReducObject(embeddings=as.matrix(combined.pca), key='PC_', assay='RNA')
  message('[scrublet atac] Finding KNN...')
  combined.seurat = combined.seurat %>%
    Seurat::L2Dim(reduction='pca') %>%
    Seurat::FindNeighbors(reduction='pca.l2', dims=dims, k=kadj, compute.SNN = FALSE) # nn.eps = 0.25    
  # From KNN, calculate doublet likelihood as defined in Scrublet paper
  message('[scrublet atac] Calculating doublet neighbors...')
  doublet_mask = ifelse(combined.seurat@meta.data$doublet, 1, 0)
  doublet_neighbors = Matrix::rowSums(safe_column_multiply(combined.seurat@graphs$RNA_nn, doublet_mask))
  message('[scrublet atac] Finalizing doublet likelihoods...')
  doublet_score = doublet_neighbors / kadj
  q = (doublet_neighbors + 1)/ (kadj + 2)
  doublet_likelihood = q * estimated_doublet_rate / fraction_sim_doublets / (1 - estimated_doublet_rate - q * (1 - estimated_doublet_rate - estimated_doublet_rate / fraction_sim_doublets))
  # Return Seurat object with doublet likelihood as an extra column
  result = data.frame(cell=rownames(combined.seurat@meta.data), 'doublet_score'=doublet_score, 'doublet_likelihood'=doublet_likelihood, 'simulated_doublet'=combined.seurat@meta.data$doublet)
  rownames(result) = result$cell
  return(result)
}
