clear
close all
clc

mkdir('./SCARDEC_results');
%% setup the key parameters.
% number of sample points for each STFs
N_pts = 100;

% normalize the STFs with the moment (integral of STF)
STF_normalization = 'integral'; % Another option using maximum: 'max_amp'

% DTW stretching constraint ratio
DTW_constraint_r = 1;

% Number of hierarchical clusters
N_clusters = 20;

% minimum prominent peak threshold for groups
MinPeakProminence_TH = 0.1; % fraction of global maximum

%% add the functions to path
addpath('./functions');

%% load and format the original results
% % This part need to be modified based on your own STF data format.
% % Can refer to the source code of read_scardec_stf.m and read_simulation.m
% % for details.
% % For testing clustering and results reproduction, can directly load my
% % processed .mat files, and run the clustering.

% mkdir('./processed_STFs')
% Scardec_Data = read_scardec_stf(['/Users/Yin9xun/Dropbox/STF_DTW_clustering/SCARDEC_original'],N_pts);
% save(['./processed_STFs/scardec_processed_stfs.mat'],'Scardec_Data');
% keyboard



%% Load the processed STF for clustering
pause(3)
disp(['   ']);
disp(['===== Loading processed SCARDEC STFs =====']);
disp(['   ']);
load('./processed_STFs/scardec_processed_stfs.mat','Scardec_Data');

%% Normalize the STFs
disp(['   ']);
disp(['===== Normalizing =====']);
disp(['   ']);
normalized_series = series_normalization(Scardec_Data.All_STFs,STF_normalization);

%% DTW clustering, stretching
disp(['   ']);
disp(['===== DTW clustering =====']);
disp(['   ']);
% FUNCTION 1: calculate the dtw distance matrix for clustering
dtw_dist = calculate_dtw_distance(normalized_series, DTW_constraint_r);

% FUNCTION 2: hierarchical clustering into n_cluster of clusters
cluster_label=hierarchical_clustering(dtw_dist, N_clusters);

% FUNCTION 3: find the center event and stretch all other stfs within each cluster
dtw_stretched_stfs = dtw_stretching(normalized_series, dtw_dist, DTW_constraint_r, cluster_label);

%% Grouping STFs based on the prominent peaks
disp(['   ']);
disp(['===== Group with prominent peaks =====']);
disp(['   ']);
group_label=prominent_peak_grouping(dtw_stretched_stfs.stretched_reference_STF,cluster_label,MinPeakProminence_TH);

%% Output clustering results for further plotting
disp(['   ']);
disp(['===== Output results =====']);
disp(['   ']);
save('./SCARDEC_results/SCARDEC_DTW_results.mat')

%% Plotting and figure reproduction
plotting_scardec_results

