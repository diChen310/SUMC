library(BayesSpace)
library(scater)
library(monocle3)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(cluster)
library(Seurat)
library(ggpubr)
library(cowplot)
library(reshape2)

K <- 100
alpha <- 0.5

set.seed(1222)
marker_sim <- function(st.i,cluster.res,K=100){
  expr.d = st.i@assays$Spatial@counts
  
  rowData <- data.frame(gene_id = rownames(expr.d),gene_short_name = rownames(expr.d),
                        sumC=Matrix::rowSums(expr.d),
                        noC = Matrix::rowSums(expr.d!=0))
  row_big <- rowData[rowData$sumC>100 & rowData$noC > 10,'gene_id']                            
  rownames(rowData) = rowData$gene_id     
  
  
  ####change into monocle3 object
  cds <- new_cell_data_set(as(expr.d[row_big,], "sparseMatrix"),
                           cell_metadata = st.i@meta.data,
                           gene_metadata = rowData[row_big,])
  cluster.res = cluster.res[Matrix::colSums(exprs(cds)) != 0]
  cds <- cds[,Matrix::colSums(exprs(cds)) != 0]
  cds <- estimate_size_factors(cds)
  cds[['bs_cluster']] = cluster.res
  cds = cds[!grepl('^MT|^RP',rownames(cds)),]
  markers_i <- top_markers(cds,group_cells_by = 'bs_cluster',
                           genes_to_test_per_group= K, cores=32)
  
  return(markers_i)
  
}

####st.objs is a list of seurat objects

####Try different q in each sample####


for( i in 1:length(st.objs)){
  
  ####read data
  st.i = st.objs[[i]]
  
  st.i = st.i[,st.i$nFeature_Spatial>0]
  image.i = st.i@images[[1]]
  image.i = image.i@coordinates
  st.i = AddMetaData(st.i,image.i)
  
  ####change into sce
  counts <- Matrix::as.matrix(st.i@assays$Spatial@counts)
  
  rowData <- data.frame(gene_id = rownames(counts),gene_name = rownames(counts),
                        sumC=apply(counts,1,sum))
  row_big <- rowData[rowData$sumC>5,'gene_id']
  
  colData <- st.i@meta.data
  colData$spot <- rownames(colData)  
  colData$in_tissue = colData$tissue
  
  sce <- SingleCellExperiment(assays = list(counts = counts[row_big,]), 
                              rowData = rowData[row_big,], colData = colData)
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  
  #####BayesSpace analysis#####
  sce <- spatialPreprocess(sce,  n.PCs = 15)
  
  #####Try K####
  sce <- spatialCluster(sce, q = 10,platform="Visium", d=10,nrep=3000,init.method="mclust", model="t")
  st.i$bs_perS_k10 = paste0('C',sce$spatial.cluster)
  
  sce <- spatialCluster(sce, q = 15,platform="Visium", d=10,nrep=3000,init.method="mclust", model="t")
  st.i$bs_perS_k15 = paste0('C',sce$spatial.cluster)
  
  sce <- spatialCluster(sce, q = 20, platform="Visium", d=10,nrep=3000,init.method="mclust", model="t")
  st.i$bs_perS_k20 = paste0('C',sce$spatial.cluster)
  
  sce <- spatialCluster(sce, q = 30,  platform="Visium", d=10,nrep=3000,init.method="mclust", model="t")
  st.i$bs_perS_k30 = paste0('C',sce$spatial.cluster)
  
  sce <- spatialCluster(sce, q = 40, platform="Visium", d=10,nrep=3000,init.method="mclust", model="t")
  st.i$bs_perS_k40 = paste0('C',sce$spatial.cluster)
  
  sce <- spatialCluster(sce, q = 50, platform="Visium", d=10,nrep=3000,init.method="mclust", model="t")
  st.i$bs_perS_k50 = paste0('C',sce$spatial.cluster)
  
  
  #write.csv(st.i@meta.data,file = paste0('./variables/bayesSpace_perS/bs_res_SS',i,'.csv'))
  
  
  print(paste('Finished for',i,'!!!!!!!!!!!!!!!'))
  
  
  
  
  sim_res_k10 <- marker_sim(st.i,st.i$bs_perS_k10)
  
  marker_k10 = sim_res_k10 %>% filter(marker_test_q_value < 0.01 & fraction_expressing > 0.3)
  
  sim_res_k15 <- marker_sim(st.i,st.i$bs_perS_k15)
  
  marker_k15 = sim_res_k15 %>% filter(marker_test_q_value < 0.01 & fraction_expressing > 0.3)
  
  sim_res_k20 <- marker_sim(st.i,st.i$bs_perS_k20)
  
  marker_k20 = sim_res_k20 %>% filter(marker_test_q_value < 0.01 & fraction_expressing > 0.3)
  
  sim_res_k30 <- marker_sim(st.i,st.i$bs_perS_k30)
  
  marker_k30 = sim_res_k30 %>% filter(marker_test_q_value < 0.01 & fraction_expressing > 0.3)
  
  sim_res_k40 <- marker_sim(st.i,st.i$bs_perS_k40)
  
  marker_k40 = sim_res_k40 %>% filter(marker_test_q_value < 0.01 & fraction_expressing > 0.3)
  
  
  sim_res_k50 <- marker_sim(st.i,st.i$bs_perS_k50)
  
  marker_k50 = sim_res_k50 %>% filter(marker_test_q_value < 0.01 & fraction_expressing > 0.3)
  
  # 
  marker_unique_no <- c(length(unique(marker_k10$gene_id)),
                        length(unique(marker_k15$gene_id)),
                        length(unique(marker_k20$gene_id)),
                        length(unique(marker_k30$gene_id)),
                        length(unique(marker_k40$gene_id)),
                        length(unique(marker_k50$gene_id))
  )
  
  
  
  print(marker_unique_no)
  
  markers_try <- list(sim_res_k10,
                      sim_res_k15,
                      sim_res_k20,
                      sim_res_k30,
                      sim_res_k40,
                      sim_res_k50
  )
  
  names(markers_try) <- cluster_tags
  sel_cluster = cluster_tags[marker_unique_no == max(marker_unique_no)]
  
  if(length(sel_cluster)>1){
    sel_cluster <- sel_cluster[1]
  }
  print(sel_cluster)
  st.i$best_cluster = st.i@meta.data[,sel_cluster]
  
  marker_best_i <- markers_try[[sel_cluster]]
  print(setdiff(unique(marker_best_i$cell_group),unique(st.i$best_cluster)))
  
  write.csv(st.i@meta.data,file = paste0('./bs_res_SSS',i,'.csv'))
  write.csv(marker_best_i,file = paste0('./bs_res_marker_top_',i,'.csv'))
  
  print(paste('Item',i,paste(marker_unique_no,collapse = ' ')))
  
}


####Merge individual local clusters####


markers_all = data.frame()
for(i in 1:length(st.objs)){
  
  markers_i = read.csv(file = paste0('./bs_res_marker_top_',i,'.csv'))
  
  print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
  print(i)
  print(length(unique(markers_i$cell_group)))
  markers_i$SampleID = unique(as.character(st.objs[[i]]$orig.ident))
  markers_all = rbind(markers_all,markers_i)
}
markers_all$bs_cluster = paste(markers_all$SampleID,'_',markers_all$cell_group)

markers_cds_all_sig <- markers_all[markers_all$marker_test_q_value<0.01,]

markers_cds_all_sig_top50 <- group_by(markers_cds_all_sig,bs_cluster) %>% top_n(K*alpha,specificity)

markers_cds_all_sig_top50 <- as.data.frame(markers_cds_all_sig_top50)

clusters_all <- unique(markers_cds_all_sig$bs_cluster)

####Matrix of gene similarity####
matrix_sim <- matrix(0,nrow = length(clusters_all),ncol = length(clusters_all),
                     dimnames = list(clusters_all,clusters_all))

for(i in 1:(length(clusters_all)-1)){
  cluster_i = clusters_all[i]
  markers_sig_i = markers_cds_all_sig_top50[markers_cds_all_sig_top50$bs_cluster == cluster_i,'gene_id']
  
  for(j in (i+1):length(clusters_all)){
    cluster_j = clusters_all[j]
    markers_sig_j = markers_cds_all_sig_top50[markers_cds_all_sig_top50$bs_cluster == cluster_j,'gene_id']
    
    inter_ij = length(intersect(markers_sig_i,markers_sig_j))
    union_ij = length(union(markers_sig_i,markers_sig_j))
    
    sim_ij = inter_ij/union_ij
    
    matrix_sim[i,j]=sim_ij
    matrix_sim[j,i] = sim_ij
  }
  print(paste(
    'Finished for',
    cluster_i
  ))
}


####gene signature####

genes_sig <- unique(markers_cds_all_sig$gene_id)

markers_cds_all_sig_temp <- markers_cds_all_sig[markers_cds_all_sig$gene_id %in% genes_sig,]
mat_markers_temp <- acast(markers_cds_all_sig_temp,bs_cluster ~ gene_id , value.var = 'marker_score',fill=0)
mat_markers_temp <- mat_markers_temp[clusters_all,]
mat_markers_temp[mat_markers_temp > 0] = 1

all.equal(rownames(matrix_sim),rownames(mat_markers_temp))

genes_sum <- apply(as.matrix(mat_markers_temp), 2, sum)

genes_l <- names(genes_sum)[genes_sum>2]

mat_gene <- mat_markers_temp[,genes_l]
####output the files for input of SUMC
write.csv(mat_gene,file='./bestk_try10_input_matrix_0.csv')
write.csv(matrix_sim,file='./bestk_try10_adj_matrix_0.csv')
