## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, 
  fig.height=5
)

## ----setup--------------------------------------------------------------------
library(SCOPRO)
required <- c("CIARA")
if (!all(unlist(lapply(required, function(pkg) requireNamespace(pkg, quietly = TRUE)))))
  knitr::opts_chunk$set(eval = FALSE)

## ---- eval = FALSE------------------------------------------------------------
#  current_wd <- getwd()
#  url = "https://hmgubox2.helmholtz-muenchen.de/index.php/s/EHQSnjMJxkR7QYT/download/SCOPRO.zip"
#  destfile <- paste0(current_wd,"/SCOPRO.zip")
#  download.file(url, destfile, quiet = FALSE)
#  unzip(destfile, exdir=current_wd)

## ---- eval = FALSE------------------------------------------------------------
#  setwd(paste0(current_wd,"/SCOPRO"))
#  load(file='mayra_dati_raw_0.Rda')
#  mayra_seurat_0=cluster_analysis_integrate_rare(mayra_dati_raw_0,"Mayra_data_0",0.1,5,30)
#  norm_es_vitro=as.matrix(GetAssayData(mayra_seurat_0, slot = "data",assay="RNA"))
#  cluster_es_vitro=as.vector(mayra_seurat_0$RNA_snn_res.0.1)

## ---- eval = FALSE------------------------------------------------------------
#  setwd(paste0(current_wd,"/SCOPRO"))
#  load(file="seurat_genes_published_mouse.Rda")
#  
#  norm_vivo <- as.matrix(GetAssayData(seurat_genes_published_mouse, slot = "data",assay="RNA"))
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  
#  DefaultAssay(seurat_genes_published_mouse) <- "RNA"
#  cluster_mouse_published <- as.vector(seurat_genes_published_mouse$stim)
#  
#  
#  relevant_stages <- c("Late_2_cell", "epiblast_4.5", "epiblast_5.5", "epiblast_6.5")
#  
#  DefaultAssay(seurat_genes_published_mouse) <- "RNA"
#  
#  markers_first_ESC_small <- CIARA::markers_cluster_seurat(seurat_genes_published_mouse[,cluster_mouse_published%in%relevant_stages],cluster_mouse_published[cluster_mouse_published%in%relevant_stages],names(seurat_genes_published_mouse$RNA_snn_res.0.2)[cluster_mouse_published%in%relevant_stages],10)
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  
#  markers_mouse <- as.vector(markers_first_ESC_small[[3]])
#  stages_markers <- names(markers_first_ESC_small[[3]])
#  
#  ## Keeping only the genes in common between in vitro and in vivo datasets
#  stages_markers <- stages_markers[markers_mouse %in% row.names(norm_es_vitro)]
#  
#  markers_small <- markers_mouse[markers_mouse %in% row.names(norm_es_vitro)]
#  names(markers_small) <- stages_markers

## ----eval = FALSE-------------------------------------------------------------
#  
#  
#  marker_result <- select_top_markers(relevant_stages, cluster_mouse_published, norm_vivo, markers_small, max_number = 100, threshold = 0.1)
#  marker_all <- marker_result[[1]]
#  marker_stages <- marker_result[[2]]
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  
#  marker_stages_filter <- filter_in_vitro(norm_es_vitro,cluster_es_vitro ,marker_all, fraction = 0.10, threshold = 0)
#  
#  analysis_2cell <- SCOPRO(norm_es_vitro,norm_vivo,cluster_es_vitro,cluster_mouse_published,"Late_2_cell",marker_stages_filter, threshold = 0.1, number_link = 1, fold_change = 3, threshold_fold_change = 0.1 ,marker_stages, relevant_stages)
#  
#  
#  
#  #png("/Users/gabriele.lubatti/Downloads/SCOPRO_1.png")
#  plot_score(analysis_2cell, marker_stages, marker_stages_filter, relevant_stages, "Late_2_cell", "Final score", "Cluster", "Late_2_cell")
#  #dev.off()
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  
#  common_genes <- select_common_genes(analysis_2cell, marker_stages, relevant_stages, "Late_2_cell", cluster_es_vitro, "2")
#  no_common_genes <- select_no_common_genes(analysis_2cell, marker_stages, relevant_stages, "Late_2_cell", cluster_es_vitro, "2")
#  
#  
#  
#  
#  all_genes <- c(no_common_genes[1:4], common_genes[1:10])
#  all_genes_label <- c(paste0(no_common_genes[1:4], "-no_conserved"), paste0(common_genes[1:10], "-conserved"))
#  
#  
#  
#  
#  
#  rabbit_plot <- plot_score_genes(all_genes, "Mouse ESC", "Mouse vitro", norm_es_vitro,norm_vivo[ , cluster_mouse_published=="Late_2_cell"],cluster_es_vitro, cluster_mouse_published[cluster_mouse_published == "Late_2_cell"], all_genes_label, 7, 10, "Late_2_cell")
#  #png("/Users/gabriele.lubatti/Downloads/SCOPRO_2.png")
#  rabbit_plot
#  #dev.off()
#  
#  
#  
#  

## -----------------------------------------------------------------------------
utils::sessionInfo()

