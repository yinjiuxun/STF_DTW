function group_label=prominent_peak_grouping(stretched_reference_STF,cluster_label,MinPeakProminence_TH)
% further group the STFs based on the number of prominent peaks of the
% center event and attach the information of prominent peaks.
% group_label=prominent_peak_grouping(stretched_reference_STF,cluster_label,MinPeakProminence_TH)
% INPUTs:
% stretched_reference_STF: matrix of stretching center events
% (dtw_stretched_stfs.stretched_reference_STF) from DTW_cluster function
% cluster_label: vectors fo DTW clustering labels of all the STFs 
% MinPeakProminence_TH: prominent peak threshold (fraction of global
% maximum).
% OUTPUT:
% group_label: group label for all STFs



% number of clusters and STFs
n_cluster = length(unique(cluster_label));
n_stf = length(cluster_label);

% find the prominent peaks and label event group
group_label = zeros(n_stf,1);
for I_cluster=1:n_cluster
    % use the number of peaks of the center event as the group number
    center_STF = stretched_reference_STF(I_cluster,:);
    [STF_PKS,~]=findpeaks(center_STF,...
        'MinPeakProminence',MinPeakProminence_TH*max(center_STF));
    
    % set the group label for each STF
    group_label(cluster_label == I_cluster) = length(STF_PKS);
end
    
    
    



