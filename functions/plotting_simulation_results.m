clear
close all
clc

% Make output directory
output_dir = './Figures';
mkdir(output_dir);

%% Colormap for the groups
GROUP_CMP=[32 100 154;
    60 173 162;
    245 212 93;
    238 86 60]/255;

%% loop over each Dc
DCs=[0.05 0.1 0.2 0.4 0.8 1.6];


%% load STFs
for dc = DCs
    
    load(['./Simulation_results/simulation_DTW_Dc_' num2str(dc) '_results.mat'])
    
    N_pts = size(Simulated_Data.All_STFs,2);
    n_stf = size(Simulated_Data.All_STFs,1);
    
    % getting the source parameters
    All_STFs = Simulated_Data.All_STFs;
    Duration = Simulated_Data.Event_info.Duration;
    radiated_energy = Simulated_Data.Event_info.radiated_energy;
    total_energy = Simulated_Data.Event_info.total_energy;
    
    radiation_efficiency = radiated_energy./total_energy;
    
    MinPeakProminence_TH=0.1; % percent of maximum
    

    % reorganize the group label (cutoff complexity only to 4 (>=4 -> 4) groups)
    Complexity_number=group_label;
    Complexity_number(Complexity_number>4) = 4;
    Complexity_label = {'G1','G2','G3','G4'};
    
    %% ===================== Reproducing Figure 1 RIGHT half =================================
    f1=figure(1);
    subplot(133)
    H=histogram(categorical(cluster_label),'DisplayOrder','ascend');
    H.Orientation='horizontal';
    H.Normalization='probability';
    H.FaceColor='b';
    set(gca, 'box','off','YTickLabel',[],'YTick',[],'XScale','linear')
    xlim([0 0.4])
    grid on
    xlabel('Number of event')
    title('Distributions')
    
    [Num_of_event,Label_event]=histcounts(categorical(cluster_label));
    
    % Resort the label of event
    [Num_of_event,Indx_sort]=sort(Num_of_event,'descend');
    Label_event=Label_event(Indx_sort);
    
    %% Go through event and Plotting the stretched and unstretched results
    stretched_reference_STF = dtw_stretched_stfs.stretched_reference_STF;
    stretched_STF = dtw_stretched_stfs.stretched_STF;
    reference_STF_indice = dtw_stretched_stfs.reference_STF_indice;

    for ic=1:length(Label_event)
        I_cluster=str2double(cell2mat(Label_event(ic))); % value of this cluster label
        % line color index
        
        event_indx_cluster=find(cluster_label==I_cluster); % event index within each cluster
                
        [STF_PKS,STF_PKS_LOCS]=findpeaks(stretched_reference_STF(I_cluster,:),...
            'MinPeakProminence',MinPeakProminence_TH*max(stretched_reference_STF(I_cluster,:)));
        
        GROUP_NUMBER = Complexity_number(reference_STF_indice(I_cluster));
        GROUP_LABEL = Complexity_label{GROUP_NUMBER};
        GROUP_color = GROUP_CMP(GROUP_NUMBER,:);
        
        
        %Plotting distributions and STFs
        MAX_temp1=max(max(stretched_STF(event_indx_cluster,:)));
        subplot(1,12,[5 7])
        plot(1:N_pts,stretched_STF(event_indx_cluster,:)'/MAX_temp1+1-ic,'-','LineWidth',1,'Color',[0.5 0.5 0.5 0.7]);
        hold on
        plot(1:N_pts,stretched_reference_STF(I_cluster,:)/MAX_temp1+1-ic,'-k','LineWidth',2.5);
        plot(STF_PKS_LOCS,STF_PKS/MAX_temp1+1-ic,'o','MarkerSize',8,'Color',GROUP_color,'MarkerFaceColor',GROUP_color);
        
        title('Stretched STFs')
        set(gca, 'box','off','XColor','w','YColor','w','TickDir','out')
        adding_text(N_pts,ic,GROUP_color,GROUP_LABEL,event_indx_cluster)
        ylim([-19 1])
           
        % plot unstretched STFs
        subplot(1,4,1)
        Matrix_stf = normalized_series;
        MAX_temp2=max(max(Matrix_stf(event_indx_cluster,:)));
        plot((1:N_pts),Matrix_stf(event_indx_cluster,:)'/MAX_temp2+1-ic,'-','LineWidth',1,'Color',[0.5 0.5 0.5 0.7]);
        hold on
        plot((1:N_pts),stretched_reference_STF(I_cluster,:)/MAX_temp1+1-ic,'-k','LineWidth',2.5);
        adding_text(N_pts,ic,GROUP_color,GROUP_LABEL,event_indx_cluster)
        titleTEXT='Unstretched STFs';
        
        ylim([-19 1])
        title(titleTEXT);
        set(gca, 'box','off','XColor','w','YColor','w','TickDir','out')
        
        
        
    end
    
    f1.Position=[50 50 800 927];
    f1.PaperSize=f1.Position(3:4);
    
    print('-dpdf','-painters',...
        [output_dir '/Fig1_right_Simulation_STFs_Dc_' num2str(dc) '.pdf'])
    %keyboard
    pause(3)
    close(f1)
end

%% ===================== Reproducing Figure 5 =================================
%% Further grouping for each individual Dc subset, bar plot
clear group_label
complexity_cutoff = 4;
%load grouping results of SCARDEC 
load('./SCARDEC_results/SCARDEC_DTW_results.mat','group_label');

Complexity_label = {'Group 1','Group 2','Group 3','Group 4'};
group_label(group_label>=complexity_cutoff) =complexity_cutoff;
Complexity_SCARDEC = Complexity_label(group_label)';

f5 = figure(5);


DC_list = [];
GROUP_complexity = [];
% adding the SCARDEC (specified as DC_list = 999
GROUP_complexity = [GROUP_complexity; categorical(Complexity_SCARDEC)];
DC_list = [DC_list; categorical(999*ones(size(Complexity_SCARDEC)))];

for dc=DCs
    load(['./Simulation_results/simulation_DTW_Dc_' num2str(dc) '_results.mat'],'group_label')
    group_label(group_label>=complexity_cutoff) =complexity_cutoff;
    Complexity_simulation = Complexity_label(group_label)';
    DC_list = [DC_list; categorical(ones(length(Complexity_simulation),1),1,{num2str(dc)})];
    GROUP_complexity = [GROUP_complexity; categorical(Complexity_simulation)];
end

DC_list_double = double(DC_list);
GROUP_complexity_double = double(GROUP_complexity);

% caluculate all frequencies:
data = accumarray([DC_list_double GROUP_complexity_double],1);
data = data./sum(data,2);
% get the categories names:
DC_list_unique = categories(DC_list); DC_list_unique{1}='SCARDEC';
GROUP_complexity_unique = categories(GROUP_complexity);
% plotting:
bb=bar(data,0.4,'stacked','FaceColor','flat','LineWidth',2);
for ii=1:complexity_cutoff
    bb(ii).CData = ii;
end
ax = gca;
ax.XTickLabel = DC_list_unique; % set the x-axis ticks to the race names
lgd=legend(GROUP_complexity_unique,'FontName','Times','FontSize',20); % add a legend for the colors
lgd.Position = [0.8530 0.1258 0.0950 0.1325];
colormap(GROUP_CMP(1:complexity_cutoff,:)) % use nicer colors (close to your example)
ylabel('Group proportion')% set the y-axis lable
ylim([0 1.1])
xlabel('Values of simulation Dc (m)')
% some other minor fixes:
box off
ax.YGrid = 'on';

f5.Position = [68 50 1000 600];
f5.PaperSize = f5.Position(3:4);

print('-dpdf',[output_dir '/Fig5_SCARDEC_vs_Simulation_group_population.pdf']);



function adding_text(N_pt,ic,GROUP_color,GROUP_LABEL,event_indx_cluster)
text(0.99*(N_pt),1-ic+0.5,['\color[rgb]{' num2str(GROUP_color) '}' GROUP_LABEL '(' num2str(length(event_indx_cluster)) ')'],...
    'FontSize',18,'HorizontalAlignment','left','Interpreter','Tex')
end




