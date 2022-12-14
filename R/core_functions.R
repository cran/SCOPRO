

#' SCOPRO
#'
#' @param norm_vitro Norm count matrix (n_genes X n_cells) for in vitro dataset
#' @param norm_vivo Norm count matrix (n_genes X n_cells) for in vivo dataset
#' @param cluster_vitro cluster for in vitro dataset
#' @param cluster_vivo  cluster for in vivo dataset
#' @param name_vivo  name of the in vivo stage on which SCOPRO is run
#' @param marker_stages_filter  output from the function \emph{filter_in_vitro}
#' @param threshold Numeric value. For a given gene, the jaccard index between the links from the in vivo and in vitro datasets is computed. If the jaccard index is above \emph{threshold}, then the gene is considered to be conserved between the two datasets.
#' @param number_link Numeric value. For a given gene in the in vivo dataset with links above \emph{number_link}, the jaccard index between the links from in vitro and in vivo dataset is computed.
#' @param threshold_fold_change Numeric value. Above \emph{threshold} the fold change between genes is computed. Below \emph{threshold} the difference between genes is computed.
#' @param fold_change Numeric value. For a given gene, the fold change between all the other genes is computed. If fold change is above \emph{fold_change}, then there is a link with weight 1 between the two genes.
#' @param marker_stages Second element of the list given as output by the function \emph{select_top_markers}
#' @param selected_stages In vivo stages for which the markers where computed with the function \emph{select_top_markers}
#' @description The mean expression profile of \emph{marker_stages_filter} genes is computed for each cluster in the in vivo and in vitro dataset.
#' For a given cluster, a connectivity matrix is computed with number of rows and number of columns equal to the length of \emph{marker_stages_filter}. Each entry (i,j)  in the matrix can be 1 if the fold_change between gene i and gene j is above \emph{fold_change}. Otherwise is 0.
#' Finally the connectivity matrix of a given \emph{name_vivo} stage and all the clusters in the in vitro dataset are compared.
#' A gene i is considered to be conserved between \emph{name_vivo} and an in vitro cluster if the jaccard index of the links of gene i is above \emph{threshold}.
#' @return A list with five elements:
#'
#' \item{common_link}{Vector with the names of the genes conserved between \emph{name_vivo} and all the clusters in the vitro dataset}
#' \item{no_common_link}{Vector with the names of the genes not conserved between \emph{name_vivo} and  the clusters in the vitro dataset}
#' \item{link_kept}{List with the names of the genes conserved between \emph{name_vivo} and each single cluster in the vitro dataset}
#' \item{link_no_kept}{List with the names of the genes not conserved between \emph{name_vivo} and each single cluster in the vitro dataset}
#' \item{final_score}{Numeric value, given by the fraction of conserved markers of \emph{name_vivo} and each single cluster in the in vitro dataset}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @examples
#' load(system.file("extdata", "norm_es_vitro_small.Rda", package = "SCOPRO"))
#' n_es= norm_es_vitro_small
#' load(system.file("extdata", "norm_vivo_small.Rda", package = "SCOPRO"))
#' n_v = norm_vivo_small
#' load(system.file("extdata", "cluster_es_vitro_small.Rda", package = "SCOPRO"))
#' c_es=cluster_es_vitro_small
#' load(system.file("extdata", "cluster_vivo_small.Rda", package = "SCOPRO"))
#' c_v=cluster_vivo_small
#' load(system.file("extdata", "marker_stages_filter.Rda", package = "SCOPRO"))
#' m_s_f = marker_stages_filter
#' load(system.file("extdata", "marker_stages.Rda", package = "SCOPRO"))
#' m_s = marker_stages
#' stages = c("Late_2_cell","epiblast_4.5","epiblast_5.5","epiblast_6.5")
#' output_SCOPRO = SCOPRO(n_es,n_v,c_es,c_v,"Late_2_cell",m_s_f,0.1,1,3,0.1,m_s,stages)
#' plot_score(output_SCOPRO,m_s,m_s_f,stages,"Late_2_cell","Score","Cluster","2-cells")
#'
#'
#' @export SCOPRO
#'

SCOPRO <- function(norm_vitro, norm_vivo, cluster_vitro, cluster_vivo, name_vivo, marker_stages_filter, threshold = 0.1, number_link = 1, fold_change = 3, threshold_fold_change = 0.1,  marker_stages, selected_stages){
  if (sum(selected_stages %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector selected_stages")

  }

  if (sum(cluster_vivo %in% name_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector cluster_vivo")

  }

  if (sum(selected_stages %in% cluster_vivo) == 0) {
    stop("name_vivo must be one the stages present in the vector cluster_vivo")

  }

  if (length(marker_stages_filter) == 0) {
    stop("Vector marker_stages_filter has 0 length. Please provide a non-zero length vector.")

  }

  mean_2_cell <- mean_stage(norm_vivo, cluster_vivo, name_vivo, marker_stages_filter)
  connettivity_vivo <- find_connettivity(mean_2_cell, fold_change, threshold_fold_change)
  name_cluster <- levels(factor(cluster_vitro))
  link_kept <- rep(list(0), length(name_cluster))
  link_no_kept <- rep(list(0), length(name_cluster))
  final_score <- rep(list(0), length(name_cluster))
  for ( i in 1:length(name_cluster)){
    mean_cluster <- mean_stage(norm_vitro,cluster_vitro,name_cluster[i], marker_stages_filter)
    connettivity_vitro <- find_connettivity(mean_cluster,fold_change, threshold_fold_change)
    connettivity_vitro <- connettivity_vitro[row.names(connettivity_vivo),colnames(connettivity_vivo)]
    index_specific <- which(selected_stages %in% name_vivo)
    marker_specific <- marker_stages[[index_specific]]
    if (length(marker_specific) == 0){
      stop(paste0("There are no markers for the stage: ",name_vivo))
    }
    result_comparison <- comparison_vivo_vitro(connettivity_vivo, connettivity_vitro, threshold,number_link, marker_specific)
    link_kept[[i]] <- names(result_comparison[[1]])[result_comparison[[2]]==1]
    if (length(link_kept[[i]]) == 0){
      warning (paste0("There are not conserved genes between ",name_vivo, " and ",name_cluster[i]))
    }
    link_no_kept[[i]] <- names(result_comparison[[1]])[result_comparison[[2]]!=1]
    names(link_kept[[i]]) <- rep(name_cluster[i],length(link_kept[[i]]))
    names(link_no_kept[[i]]) <- rep(name_cluster[i],length(link_no_kept[[i]]))
    final_score[[i]] <- result_comparison[[3]]
    names(final_score[[i]]) <- name_cluster[i]
  }
  common_link <- Reduce(intersect, link_kept)
  if (length(common_link) == 0){
    warning (paste0("There are not conserved genes between all the in vitro cluster and ",name_vivo))
  }
  no_common_link <- Reduce(intersect, link_no_kept)
  return(list(common_link, no_common_link, link_kept, link_no_kept, final_score))
}

#' mean_stage
#' @noRd
mean_stage <- function(norm_counts, stage, name_stage, marker_stages_filter){
  norm_counts <- norm_counts[,stage == name_stage]
  mean_stage <- apply(norm_counts[marker_stages_filter, ], 1, mean)
  mean_stage <- sort(mean_stage, decreasing = T)
  return(mean_stage)
}

#' find_connettivity
#' @noRd
find_connettivity <- function(mean_stage, fold_change, threshold_fold_change){
  genes_fold_change <- rep(list(0), length(mean_stage))
  for (i in 1:length(mean_stage)){
    if (mean_stage[i] <= 0.1){
      ratio <- mean_stage[i]-mean_stage
      genes_fold_change[[i]][ratio == threshold_fold_change] <- 1
      genes_fold_change[[i]][ratio < threshold_fold_change] <- 0

    }
    else{
      ratio <- mean_stage[i] / mean_stage
      genes_fold_change[[i]][ratio <= fold_change] <- 0
      genes_fold_change[[i]][ratio > fold_change] <- 1
    }
  }

  connettivity <- matrix(unlist(genes_fold_change), byrow=TRUE, nrow=length(genes_fold_change) )
  row.names(connettivity) <- names(mean_stage)
  colnames(connettivity) <- names(mean_stage)
  return(connettivity)}


#' comparison_vivo_vitro
#' @noRd
comparison_vivo_vitro <- function(connettivity_vivo, connettivity_vitro, threshold, number_link, marker_specific){
  frac <- rep(0,length(row.names(connettivity_vivo)))
  for ( i in 1:length(row.names(connettivity_vivo))){
    index_vivo <- which(connettivity_vivo[i,]==1)
    index_vitro <- which(connettivity_vitro[i,]==1)
    if (length(index_vivo) > number_link){
      frac[i] <- length(intersect(index_vitro, index_vivo))/length(union(index_vitro, index_vivo))
    }
    else{frac[i]="Excluded"}
  }

  names(frac) <- row.names(connettivity_vivo)
  frac_1 <- frac[frac != "Excluded"]
  frac_1 <- sort(frac_1,decreasing = TRUE)
  frac_final <- frac_1
  frac_final[frac_1 > threshold] <- 1
  frac_final[frac_1 <= threshold] <- 0
  frac_specific <- sum(frac_1[names(frac_1)%in%marker_specific]>threshold)/length(frac_1[names(frac_1)%in%marker_specific])
  return(list(frac_1,frac_final,frac_specific))
}


