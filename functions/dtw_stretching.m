%% FUNCTION 3: Find the center event and stretch all other stfs within each cluster
function dtw_stretched_stfs = dtw_stretching(All_STFs, dtw_dist, DTW_constraint, cluster_label)
%% stretching the STFs based on the clustering results (center events)
% dtw_stretched_stfs = dtw_stretching(All_STFs, dtw_dist, DTW_constraint, cluster_label)
% INPUTs:
% All_STFs: STF matrix (n_stf,n_pts)
% dtw_dist: Calcuated DTW distance matrix (n_stf,n_stf)
% DTW_constraint: stretching constraint ratio for dtw, if traversal
% searching, set DTW_constraint = 1. For more details see dtw
% cluster_label: Cluster label vector for all STFs, varying from 1 to
% n_cluster
%
% OUTPUTs:
% dtw_stretched_stfs
%   is a structure with fields:
%     stretched_reference_STF: [n_cluster×n_pts double]
%       --> center event of each cluster 
%        reference_STF_indice: [n_cluster×1 double]
%       --> index of each center event
%               stretched_STF: [n_cluster×n_pts double]
%       --> Matrix with all stretched STFs

dtw_stretched_stfs = struct;
% number of STFs and number of points in the STF time series
n_stf = size(All_STFs,1);
n_pts = size(All_STFs,2);
% number of clusters
n_cluster = max(cluster_label);

% to record the reference/center STFs
stretched_reference_STF=zeros(n_cluster,n_pts);
reference_STF_indice = zeros(n_cluster,1); % index of center event in the orginal matrix
% to record the stretched STFs
stretched_STF=zeros(n_stf,n_pts);

for i=1:n_cluster
    I_cluster=find(cluster_label==i);
    
    %% strech each STF in one cluster based on the center of cluster
    if length(I_cluster)==1
        subXX_r=I_cluster;
    else
        %The center event is defined as the one with minimum median
        %distance to other event in one cluster
        TEMP_DIST=dtw_dist(I_cluster,I_cluster);
        median_dist=median(TEMP_DIST,1);
        [~,INX_min]=min(median_dist);
        subXX_r=I_cluster(INX_min); % the "center" of each cluster
    end
    
    % save the center event STF and index
    reference_STF_indice(i) = subXX_r;
    stretched_reference_STF(i,:) = All_STFs(subXX_r,:);
    
    % initialize the stretching indice for every STF in each cluster
    stretched_Iy_STF_all=zeros(length(I_cluster),n_pts);
    
    for ic=1:length(I_cluster)
        %%%%%%%%%%% stretching STF
        [~,Iy,Ix]=dtw(All_STFs(I_cluster(ic),:),All_STFs(subXX_r,:),round(DTW_constraint*n_pts));
        
        [~,IA1,~] = unique(Ix);  % C1=Ix(IA)
        
        temp_Iy=All_STFs(I_cluster(ic),Iy);
        temp_Iy=temp_Iy(IA1);
        stretched_Iy_STF_all(ic,:)=temp_Iy; % used for following plot purpose
        stretched_STF(I_cluster(ic),:)=temp_Iy; % store the stretched STF
        
    end
end

% assembling into one output structure.
dtw_stretched_stfs.stretched_reference_STF = stretched_reference_STF;
dtw_stretched_stfs.reference_STF_indice = reference_STF_indice; 
dtw_stretched_stfs.stretched_STF=stretched_STF;

end