% constraint in the dtw function (maximum step of matching)
DTW_constraint=500;

% do either maxima scaling or area scaling
Scaling_Index=2; % 0: maxima,  1: area, 2: area, stf only
if Scaling_Index==0
    txt_tail=['_DTW_' num2str(DTW_constraint) '_amp_scaling'];
elseif Scaling_Index==1
    txt_tail=['_DTW_' num2str(DTW_constraint) '_area_scaling_median'];
elseif Scaling_Index==2
    txt_tail=['_DTW_' num2str(DTW_constraint) '_area_scaling_median_STF_only'];
end


%% load STFs
load('All_stfs_500.mat','All_STFs','All_REs','All_headers');

All_STFs(isnan(All_STFs))=0;
All_REs(isnan(All_REs))=0;

% % for testing purpose
% All_STFs=All_STFs(1:100,:);
% All_REs=All_REs(1:100,:);

FILE_name=['./SCARDEC' txt_tail];
mkdir(FILE_name);


% normalize the STFs
N_stf=size(All_STFs,1);
N_pt=size(All_STFs,2);

Matrix_stf=zeros(N_stf,N_pt);
Matrix_re=zeros(N_stf,N_pt);

% extract some source information
Mg=All_headers(:,11); % magnitude
Depth=All_headers(:,9); % Depth


for istf=1:N_stf
    
    temp_STF=All_STFs(istf,:);
    temp_SAF=All_REs(istf,:);
  
    % normalized STF and SAF
    temp_STF_normal=temp_STF;
    temp_SAF_normal=temp_SAF;

    % scaling with maximum amplitude or area
    if Scaling_Index==0
        % amplitude scaling
        temp_STF_normal=temp_STF_normal/max(temp_STF_normal);
        temp_SAF_normal=temp_SAF_normal/max(temp_SAF_normal);
        
    elseif Scaling_Index==1 || Scaling_Index==2
        % area scaling
        temp_STF_normal=temp_STF_normal/(trapz(temp_STF_normal));
        temp_SAF_normal=temp_SAF_normal/(trapz(temp_SAF_normal));
    end
    
    Matrix_stf(istf,:)=temp_STF_normal;
    Matrix_re(istf,:)=temp_SAF_normal;
end

%% DTW distance between STFs and REs
DTW_STF_dist=zeros(N_stf,N_stf);
disp(['==============Clustering STFs=================']);
for II=1:N_stf
    
    if mod(II,20)==0
        disp(['------------Clustering STFs: ' num2str(II/N_stf*100) '% --------------']);
    end
    
    for JJ=1:N_stf
        [DTW_STF_dist(II,JJ),~,~]=dtw(Matrix_stf(II,:),Matrix_stf(JJ,:),DTW_constraint);
        
    end
end


DTW_RE_dist=DTW_STF_dist; 
if Scaling_Index ~=2
    disp(['==============Clustering REs===================']);
    for II=1:N_stf
        
        if mod(II,10)==0
            disp(['------------Clustering REs: ' num2str(II/N_stf*100) '% --------------']);
        end
        
        for JJ=1:N_stf
            [DTW_RE_dist(II,JJ),~,~]=dtw(Matrix_re(II,:),Matrix_re(JJ,:),DTW_constraint);
            
        end
    end
end


% use the geometrical average distance as the DTW distance
DTW_dist=sqrt(DTW_STF_dist.*DTW_RE_dist);
save([FILE_name '/DTW_distance' txt_tail '.mat'],'Scaling_Index','DTW_dist','DTW_STF_dist','DTW_RE_dist');

%% Unsupervised clustering based on the distances
% transform to pdist list
Z=squareform(DTW_dist);
Ztree = linkage(Z,'complete');
n_cluster=20;
cluster_label = cluster(Ztree,'maxclust',n_cluster);

% find the threshold distance
for i=1:n_cluster
    I_cluster=find(cluster_label==i);
    TEMP_DIST=DTW_dist(I_cluster,I_cluster);
    max_dist(i)=max(TEMP_DIST(:));
end
max_dist_all=max(max_dist(:));


f7=figure(7);
subplot(2,2,1)
dendrogram(Ztree,0,'ColorThreshold',max_dist_all);
xlim([-50 N_stf+50])
set(gca,'XTick',[])
box on
title('Dendrogram of Clustering')
ylabel('DTW distance')

figure(7)
subplot(2,2,3)
Hc=histogram(cluster_label,'BinEdges',0.5:1:n_cluster+1);
Hc.FaceColor=[0.7 0.7 0.7];

cluster_num=Hc.Values;

xlabel('Cluster #')
set(gca,'XTick',1:1:n_cluster)
title('Histogram of clusters')

%%

CMP=ones(64,3)*0.6;
f7=figure(7);
f7.Visible='off';

% to record the reference/center STFs
stretched_reference_STF=zeros(n_cluster,N_pt);
stretched_reference_RE=zeros(n_cluster,N_pt);
% to record the stretched STFs and SAFs
stretched_STF=zeros(N_stf,N_pt);
stretched_RE=zeros(N_stf,N_pt);

for i=1:n_cluster
    I_cluster=find(cluster_label==i);
    
    %% strech each STF in one cluster based on the center of cluster
    if length(I_cluster)==1
        subXX_r=I_cluster;
    else
        
        %The center event is defined as the one with minimum mean
        %distance to other event in one cluster
        TEMP_DIST=DTW_dist(I_cluster,I_cluster);
        mean_dist=mean(TEMP_DIST,1);
        median_dist=median(TEMP_DIST,1);
        
        [~,INX_min]=min(median_dist);
        subXX_r=I_cluster(INX_min); % the "center" of each cluster
    end
    
    
    stretched_Iy_STF_all=zeros(length(I_cluster),N_pt);
    stretched_Iy_SAF_all=zeros(length(I_cluster),N_pt);
    
    for ic=1:length(I_cluster)
        %%%%%%%%%%% stretching STF
        [~,Iy,Ix]=dtw(Matrix_stf(I_cluster(ic),:),Matrix_stf(subXX_r,:),DTW_constraint);
        
        [~,IA1,~] = unique(Ix);  % C1=Ix(IA)
        
        temp_Iy=Matrix_stf(I_cluster(ic),Iy);
        temp_Iy=temp_Iy(IA1);
        stretched_Iy_STF_all(ic,:)=temp_Iy; % used for following plot purpose
        stretched_STF(I_cluster(ic),:)=temp_Iy; % store the stretched STF
        
        if Scaling_Index==2 %DTW only on STF, so stretch the RE accordingly with STF index
            temp_Iy=Matrix_re(I_cluster(ic),Iy);
            temp_Iy=temp_Iy(IA1);
            stretched_Iy_SAF_all(ic,:)=temp_Iy;
            stretched_RE(I_cluster(ic),:)=temp_Iy; % store the stretched RE
        else
            %%%%%%%%%% stretching SAF seperately for STF x RE
            [~,Iy,Ix]=dtw(Matrix_re(I_cluster(ic),:),Matrix_re(subXX_r,:),DTW_constraint);
            
            [~,IA1,~] = unique(Ix);  % C1=Ix(IA)
            
            temp_Iy=Matrix_re(I_cluster(ic),Iy);
            temp_Iy=temp_Iy(IA1);
            stretched_Iy_SAF_all(ic,:)=temp_Iy; % used for following plot purpose
            stretched_RE(I_cluster(ic),:)=temp_Iy; % store the stretched RE
        end
    end
    %% Plotting the stretched results
    for ic=1:length(I_cluster)
        % line color index
        Ind_color=floor(63*(Mg(I_cluster(ic))-min(Mg))/(max(Mg)-min(Mg)))+1;
        
        %figure(7)
        % stretched STF
        MAX_temp=max(stretched_Iy_STF_all(:));
        subplot(2,2,[2 4])
        
        if ic==length(I_cluster)
            subplot(2,2,[2 4])
            plot(1:N_pt,stretched_Iy_STF_all(find(I_cluster==subXX_r),:)/MAX_temp+i*1,'-k','LineWidth',2.5);
        end
        
        hold on
        plot(1:N_pt,stretched_Iy_STF_all(ic,:)/MAX_temp+i*1,'-','LineWidth',1,'Color',[CMP(Ind_color,:) 1]);
        
    end
    %text(0.9*N_pt,i+0.5,num2str(length(I_cluster)),'FontSize',20)
    
    for ic=1:length(I_cluster)
        % line color index
        Ind_color=floor(63*(Mg(I_cluster(ic))-min(Mg))/(max(Mg)-min(Mg)))+1;
        
        %figure(7)
        % stretched SAF
        MAX_temp=max(stretched_Iy_SAF_all(:));
        subplot(2,2,[2 4])
        
        if ic==length(I_cluster)
            subplot(2,2,[2 4])
            plot((1:N_pt)+N_pt,stretched_Iy_SAF_all(find(I_cluster==subXX_r),:)/MAX_temp+i*1,'-r','LineWidth',2.5);
        end
        
        hold on
        plot((1:N_pt)+N_pt,stretched_Iy_SAF_all(ic,:)/MAX_temp+i*1,'-','LineWidth',1,'Color',[CMP(Ind_color,:) 1]);
        
        
    end
    text(0.99*(N_pt+100),i+0.5,['(' num2str(length(I_cluster)) ')'],'FontSize',20,'HorizontalAlignment','right')
    
    % record the reference/center STFs
    stretched_reference_STF(i,:)=stretched_Iy_STF_all(find(I_cluster==subXX_r),:);
    stretched_reference_RE(i,:)=stretched_Iy_SAF_all(find(I_cluster==subXX_r),:);
end



subplot(2,2,[2 4])
title('DTW Stretched STFs and REs')
ylabel('Cluster #')
set(gca,'YTick',1:n_cluster,'XTick',[]);
ylim([0 n_cluster+1])
colormap(CMP)

f7.PaperUnits='points';
f7.Position=[50 50 1193 1184];
f7.PaperSize=f7.Position(3:4);
f7.Visible='on';

%keyboard
print('-dpdf','-painters',[FILE_name '/DTW_SCARDEC_events_STFxSAF_cluster' txt_tail '.pdf'])

%pause
%

DTW_labels=categorical(cluster_label);
save([FILE_name '/clustered_STFs_x_SAF_PW_fix_phase' txt_tail '_DTW_info.mat'],...
    'Matrix_stf','DTW_labels','Matrix_re','stretched_STF','stretched_RE','stretched_reference_STF','stretched_reference_RE',...
    'Mg','Depth');




