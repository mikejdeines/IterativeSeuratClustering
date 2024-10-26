initial_clustering <- function(seurat_object){
  # Finds low-resolution initial clusters and runs the first iteration of clustering.
  require(Seurat)
  require(scCustomize)
  # Find initial clusters at low resolution
  seurat_object <- FindClusters(seurat_object, resolution = 0.01, algorithm = 4, method = "igraph")
  # Join RNA Assay Layers
  seurat_object@assays$RNA <- JoinLayers(seurat_object@assays$RNA)
  # Split object based on initial cluster membership
  objs <- SplitObject(seurat_object, "seurat_clusters")
  objs <- objs[order(names(objs))]
  # Divide initial clusters into 2 clusters and test DE genes
  samples <- lapply(objs, leiden_clustering)
  samples <- samples[order(names(samples))]
  merged_seurats <- merge(samples[[1]], samples[-1])
  merged_seurats@graphs <- seurat_object@graphs
  merged_seurats@reductions <- seurat_object@reductions
  return(merged_seurats)
}
clustering_iteration <- function(seurat_object, max_score){
  require(Seurat)
  require(scCustomize)
  Idents(seurat_object) <- seurat_object$leiden_clusters
  objs <- SplitObject(seurat_object, "leiden_clusters")
  objs <- objs[order(names(objs))]
  samples <- lapply(objs, leiden_clustering, score_limit = max_score)
  samples <- samples[order(names(samples))]
  merged_seurats <- merge(samples[[1]], samples[-1])
  merged_seurats@graphs <- seurat_object@graphs
  merged_seurats@reductions <- seurat_object@reductions
  return(merged_seurats)
}
iterative_clustering <- function(seurat_object, max_iterations = 10, max_score = 150){
  seurat_object <- initial_clustering(seurat_object)
  cluster_sizes <- data.frame()
  for (i in 1:max_iterations-1){
    old_clusters <- length(levels(seurat_object$leiden_clusters))
    print(paste("Previous Cluster Count:", old_clusters, sep = " "))
    seurat_object <- clustering_iteration(seurat_object, max_score)
    seurat_object$leiden_clusters <- as.factor(seurat_object$leiden_clusters)
    num_clusters <- length(levels(seurat_object$leiden_clusters))
    if (old_clusters - num_clusters == 0) break()
    }
  return(seurat_object)
}
leiden_clustering <- function(seurat_object, num_clusters = 2, score_limit = 150){
  # Cluster using Leiden algorithim into num_clusters
  require(Seurat)
  seurat_object$starting_clusters <- Idents(seurat_object)
  initial_resolution = 1
  seurat_object <- FindClusters(seurat_object, resolution = initial_resolution, algorithm = 4, method = "igraph")
  cluster_count <- length(levels(Idents(seurat_object)))
  while (cluster_count != num_clusters){
    print(paste("Cluster count:", cluster_count))
    if (cluster_count > num_clusters){
      initial_resolution = initial_resolution/10
      seurat_object <- FindClusters(seurat_object, resolution = initial_resolution, algorithm = 4, method = "igraph")
      cluster_count <- length(levels(seurat_object$seurat_clusters))
    }
    else{
      initial_resolution = initial_resolution*5
      seurat_object <- FindClusters(seurat_object, resolution = initial_resolution, algorithm = 4, method = "igraph")
      cluster_count <- length(levels(seurat_object$seurat_clusters))
    }
  }
  print(paste("Cluster count:", cluster_count))
  seurat_object$leiden_clusters <- Idents(seurat_object)
  seurat_object$leiden_clusters <- paste(seurat_object$starting_clusters, "_", seurat_object$leiden_clusters, sep = "")
  score <- clustering_score(seurat_object)
  if (score < score_limit){
    seurat_object$leiden_clusters <- seurat_object$starting_clusters
  }
  sizes <- data.frame(table(seurat_object$leiden_clusters))
  for (i in 1:length(sizes)){
    if (sizes$Freq[[i]] < 20){
      seurat_object$leiden_clusters <- seurat_object$starting_clusters
  }
  return(seurat_object)
}
clustering_score <- function(seurat_object){
  require(Seurat)
  de.genes <- FindMarkers(seurat_object, ident.1 = 1, ident.2 = 2, logfc.threshold = 1, min.pct = 0.25, recorrect_umi = FALSE)
  de.genes <- subset(de.genes, abs(de.genes$avg_log2FC) > 2)
  score <- -log10(de.genes$p_val_adj)
  score[score > 20] <- 20
  score <- sum(score)
  return(score)
}
