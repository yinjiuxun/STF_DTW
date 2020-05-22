%% FUNCTION 1: DTW distance between STFs or REs
function dtw_dist = calculate_dtw_distance(All_STFs, DTW_constraint)
%% calculate the dtw distance matrix
% dtw_dist = calculate_dtw_distance(All_STFs, DTW_constraint)
% INPUTs:
% All_STFs: STF matrix (n_stf,n_pts)
%
% DTW_constraint: stretching constraint ratio for dtw, if traversal
% searching, set DTW_constraint = 1. For more details see dtw
%
% OUTPUTs:
% dtw_dist: Symmetric matrix of dtw distance for each STF pair (n_stf,n_stf)
%

n_stf = size(All_STFs,1);
n_pts = size(All_STFs,2);
dtw_dist=zeros(n_stf,n_stf);
for II=1:n_stf
    
    if mod(II,20)==0
        disp(['------------Clustering STFs: ' num2str(II/n_stf*100) '% --------------']);
    end
    
    for JJ=II:n_stf
        [dtw_dist(II,JJ),~,~]=dtw(All_STFs(II,:),All_STFs(JJ,:),round(DTW_constraint*n_pts));
        
    end
end
dtw_dist = dtw_dist + dtw_dist';
end