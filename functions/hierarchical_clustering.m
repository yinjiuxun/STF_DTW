%% FUNCTION 2: Hierachical clustering with single linkage
function cluster_label=hierarchical_clustering(dtw_dist, n_cluster)
%% clustering the stfs based on the dtw distance matrix 
% cluster_label=hierarchical_clustering(dtw_dist, n_cluster)
% INPUTs: 
% dtw_dist: Calcuated DTW distance matrix (n_stf,n_stf)
% n_cluster: Number of clusters for the hierarchical clustering
% OUTPUTs:
% cluster_label: Cluster label vector for all STFs, varying from 1 to
% n_cluster

% transform to pdist list
Z=squareform(dtw_dist);
% build up the linkage tree
Ztree = linkage(Z,'complete');
cluster_label = cluster(Ztree,'maxclust',n_cluster);

end