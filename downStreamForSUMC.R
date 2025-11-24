
library(ConsensusClusterPlus)
library(cluster)
library(pheatmap)
####consensus clustering based on embeddings####
clusters_all <- unique(markers_cds_all_sig$bs_cluster)
mlp_emp <- read.csv(file = './pythonOut/SUMC_output_embeddings.csv',row.names = 1)
rownames(mlp_emp) <- rownames(matrix_sim)


dist_emp <- dist(mlp_emp)


cluster_labels_k30_cc <- ConsensusClusterPlus(d=t(mlp_emp),maxK = 30,clusterAlg = 'pam')

silhouette_scores <- sapply(10:30, function(k) {
  cluster_labels <- cluster_labels_k30_cc[[k]][["consensusClass"]]
  silhouette_score <- mean(silhouette(cluster_labels, dist_emp)[, 3])
  return(silhouette_score)
})

plot(10:30, silhouette_scores, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters (k)", ylab = "silhouette coefficient",
     main = "silhouette coefficient")
	 
cluster_k19_cc <- cluster_labels_k30_cc[[19]][["consensusClass"]]


res_sc = data.frame(orig.ident = unlist(lapply(clusters_all,function(a){
  strsplit(a,' _ ')[[1]][1]
})),
cluster_k30= paste('K',cluster_k30_cc),
cluster_k19= paste('M',cluster_k19_cc),
row.names = clusters_all)


####assign to each st####

for (i in 1:length(st.objs)) {
  
  st.i = st.objs[[i]]
  meta.i = read.csv(file = paste0('./variables/bayesSpace_perS/bs_res_SSS',i,'.csv'),row.names = 1)
  meta.i$bs_cluster = paste(meta.i$orig.ident,'_',meta.i$best_cluster)
  #meta.i$cluster_k10 = res_sc[meta.i$bs_cluster,'cluster_k10']
  meta.i$cluster_k30 = res_sc[meta.i$bs_cluster,'cluster_k30']
  meta.i$cluster_k19 = res_sc[meta.i$bs_cluster,'cluster_k19']
  
  st.i = AddMetaData(st.i,meta.i)
  st.objs[[i]] = st.i
}
