clear
close all
clc

% read the SCARDEC source time function, preprocessing and extract some
% of the source parameters
read_SCARDEC_STF_all;

% apply the DTW clustering to the STFs
close all
SCARDEC_STF_RE_clusters_DTW;

% plot the clustering results and correlation with the source parameters
close all
plot_DTW_500_no_constraint_sorted_by_number

%% Output a text file with source parameters and clustering results
TEMP_MTX=[Num_complexity,double(DTW_labels),All_headers];

fileID = fopen('stf_stats_yin.dat','wt');
for ii = 1:size(TEMP_MTX,1)
    fprintf(fileID,'%g\t',TEMP_MTX(ii,:));
    fprintf(fileID,'\n');
end
fclose(fileID);

%% Highlight some of the special events using the time

TEMP_MTX=load('stf_stats_yin.dat');
% year month day hour min second
EVENT_date=[2012 4 11 8 38 36;   % 2012 Sumatra
    2016 11 13 11 2 56;  % 2016 Kaikoura
    1992 6  28 11 57 34; % 1992 Landers
    2013 9 24 11 29 47;  % 2013 Pakistan
    
    2010 2 27  6 34 11; % 2010 Chile
    2011 3  11 5 46 24; % 2011 Tohoku
    2015 4 25 6 11 25; %NEPAL
    2015 9 16 22 54 32; %NEAR COAST OF CENTRAL CHILE
    2017 9 8 4 49 19]; % NEAR COAST OF CHIAPAS MEXICO

II=zeros(1,size(EVENT_date,1));
for ii_event=1:size(EVENT_date,1)
    II(ii_event)=find(TEMP_MTX(:,3)==EVENT_date(ii_event,1) & TEMP_MTX(:,4)==EVENT_date(ii_event,2) ...
        & TEMP_MTX(:,5)==EVENT_date(ii_event,3) & TEMP_MTX(:,6)==EVENT_date(ii_event,4) ...
        & TEMP_MTX(:,7)==EVENT_date(ii_event,5) & TEMP_MTX(:,8)==EVENT_date(ii_event,6));
    
end

TEMP_MTX(II,1:2)