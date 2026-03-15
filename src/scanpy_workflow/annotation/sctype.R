#!/usr/bin/env Rscript
# sctype.R — scType-style automated cell type annotation for Nextflow DSL2
#
# Usage:
#   Rscript sctype.R \
#     --input   <clustered.h5ad> \
#     --markers <markers.json> \
#     --groupby <leiden> \
#     --output  <annotated.h5ad>
#
# markers.json format:
#   { "T cells": { "pos": ["CD3D", "CD3E"], "neg": [] }, ... }
#
# Algorithm:
#   For each cluster C and each cell type T:
#     score(C, T) = mean_expr(pos_genes in C) - mean_expr(neg_genes in C)
#   Assign highest-scoring cell type (score > 0) or "Unknown" to each cluster.
#   All cells in the same cluster receive the same annotation.
#   Result stored in colData$sctype_annotation.
#
# Requires: r-jsonlite (conda), zellkonverter (bioconda)
# Optional: remotes::install_github('IanevskiAleksandr/sc-type') for full scType

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(jsonlite)
  library(Matrix)
})

parse_args_sctype <- function() {
  option_list <- list(
    make_option("--input",   type = "character", help = "Input .h5ad path"),
    make_option("--markers", type = "character", help = "Markers JSON path"),
    make_option("--groupby", type = "character", default = "leiden",
                help = "obs column with cluster assignments [default: leiden]"),
    make_option("--output",  type = "character", help = "Output .h5ad path")
  )
  parser <- OptionParser(option_list = option_list)
  opts   <- parse_args(parser)
  if (is.null(opts$input) || is.null(opts$markers) || is.null(opts$output)) {
    print_help(parser)
    stop("--input, --markers, and --output are required", call. = FALSE)
  }
  opts
}

#' Score each cluster against each cell type from the markers list.
#' @param expr_mat  genes x cells expression matrix.
#' @param clusters  Named character vector mapping cell -> cluster label.
#' @param markers   Named list: cell_type -> list(pos = ..., neg = ...).
#' @return data.frame: cluster, cell_type, score.
score_clusters <- function(expr_mat, clusters, markers) {
  genes_in_mat <- rownames(expr_mat)
  cluster_ids  <- unique(clusters)
  records <- vector("list", length(cluster_ids) * length(markers))
  k <- 1L
  for (cl in cluster_ids) {
    cell_idx <- which(clusters == cl)
    cl_means <- rowMeans(expr_mat[, cell_idx, drop = FALSE])
    for (ct in names(markers)) {
      pos_g     <- intersect(markers[[ct]]$pos, genes_in_mat)
      neg_g     <- intersect(markers[[ct]]$neg, genes_in_mat)
      pos_score <- if (length(pos_g) > 0) mean(cl_means[pos_g]) else 0
      neg_score <- if (length(neg_g) > 0) mean(cl_means[neg_g]) else 0
      records[[k]] <- data.frame(cluster = cl, cell_type = ct,
                                 score = pos_score - neg_score,
                                 stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  do.call(rbind, records)
}

run_sctype <- function(input_path, markers_path, groupby, output_path) {
  message("sctype.R: reading ", input_path)
  sce <- readH5AD(input_path, reader = "R")

  if (!groupby %in% names(colData(sce))) {
    stop("colData does not contain column '", groupby, "'. ",
         "Available: ", paste(names(colData(sce)), collapse = ", "))
  }

  expr_mat <- if ("log_norm" %in% assayNames(sce)) {
    as.matrix(assay(sce, "log_norm"))
  } else {
    message("sctype.R: 'log_norm' assay not found; using 'X'")
    as.matrix(assay(sce, "X"))
  }

  message("sctype.R: loading markers from ", markers_path)
  # simplifyVector = TRUE (default): JSON arrays become atomic character vectors.
  # simplifyVector = FALSE would return list objects that break intersect().
  markers <- fromJSON(markers_path)
  for (ct in names(markers)) {
    # Coerce to character; handles NULL (absent key) and numeric/list edge cases
    if (is.null(markers[[ct]]$pos)) markers[[ct]]$pos <- character(0)
    if (is.null(markers[[ct]]$neg)) markers[[ct]]$neg <- character(0)
    markers[[ct]]$pos <- as.character(markers[[ct]]$pos)
    markers[[ct]]$neg <- as.character(markers[[ct]]$neg)
  }

  clusters <- setNames(as.character(colData(sce)[[groupby]]), colnames(sce))
  message("sctype.R: scoring ", length(unique(clusters)), " clusters against ",
          length(markers), " cell types")

  scores_df <- score_clusters(expr_mat, clusters, markers)

  # Per cluster: pick highest-scoring cell type (or "Unknown" if max <= 0)
  best_df <- do.call(rbind, lapply(split(scores_df, scores_df$cluster), function(df) {
    best_idx <- which.max(df$score)
    label    <- if (df$score[best_idx] > 0) df$cell_type[best_idx] else "Unknown"
    data.frame(cluster = df$cluster[1], annotation = label, stringsAsFactors = FALSE)
  }))

  cluster_to_annot      <- setNames(best_df$annotation, best_df$cluster)
  sce$sctype_annotation <- cluster_to_annot[clusters]

  message("sctype.R: annotation summary:")
  print(table(sce$sctype_annotation))

  message("sctype.R: writing ", output_path)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeH5AD(sce, output_path)
  message("sctype.R: done")
}

main <- function() {
  opts <- parse_args_sctype()
  run_sctype(opts$input, opts$markers, opts$groupby, opts$output)
}

if (!interactive()) main()
