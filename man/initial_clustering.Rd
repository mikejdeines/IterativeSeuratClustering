% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iterative_clustering.R
\name{initial_clustering}
\alias{initial_clustering}
\title{Performs an initial clustering iteration on a Seurat object.
Runs Leiden clustering at low resolution to find the initial clusters and runs one clustering iteration.}
\usage{
initial_clustering(seurat_object)
}
\arguments{
\item{seurat_object}{a normalized, integrated Seurat object}
}
\value{
a Seurat object with initial clusters in the "leiden_clusters" slot
}
\description{
Performs an initial clustering iteration on a Seurat object.
Runs Leiden clustering at low resolution to find the initial clusters and runs one clustering iteration.
}
