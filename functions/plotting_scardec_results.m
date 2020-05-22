clear
close all
clc

%% Load the results from DTW clustering
addpath('./functions');
load('./SCARDEC_results/SCARDEC_DTW_results.mat')

N_pts = size(Scardec_Data.All_STFs,2);
n_stf = size(Scardec_Data.All_STFs,1);
total_RE_all = zeros(n_stf,1);

% Make output directory
output_dir = './Figures';
mkdir(output_dir);

%% Looking at the correlations between source parameters and group (only look at 4 groups)
All_STFs = Scardec_Data.All_STFs;
Lat = Scardec_Data.Event_info.Lat;
Lon = Scardec_Data.Event_info.Lon;
Depth = Scardec_Data.Event_info.Depth;
Moment0 = Scardec_Data.Event_info.Moment; % seismic moment
Mw = Scardec_Data.Event_info.Mw; % moment magnitude
gaussian_subevent_number = Scardec_Data.Event_info.gaussian_subevent_number; % Gaussian subevent number
FocalMechanism = Scardec_Data.Event_info.FocalMechanism; % focal mechanism [strike1,dip1,rake1,strike2,dip2,rake2]


T_duration = Scardec_Data.T_duration;

fm_vector=focal_mechanism_vector(Scardec_Data.Event_info.FocalMechanism(:,[3 6]));

% calculate the integral of squared slip acceleration function 
for i_stf=1:length(Depth)
    Time_temp = linspace(0,T_duration(i_stf),N_pts);
    RE_temp = Scardec_Data.All_REs(i_stf,:);
    total_RE_all(i_stf)=trapz(Time_temp,RE_temp);  
end

%Designa is stress drop in MPa calculated based on Eshelby 1957, assuming that size (Brune 1970) 
TempA=load('prem_model.txt'); 
depth_prem=TempA(:,1);
Vs_prem=TempA(:,4);
Vp_prem=TempA(:,3);
Rho_prem=TempA(:,2);

k=0.32;
Vs=interp1(depth_prem,Vs_prem,Depth)*1e3;
Vp=interp1(depth_prem,Vp_prem,Depth)*1e3;
Rho=interp1(depth_prem,Rho_prem,Depth)*1e3;
% Scaled energy ER/M0
total_RE_all_scaled=(1/15/pi./Rho./Vp.^5 + 1/10/pi./Rho./Vs.^5).*total_RE_all./Moment0;

Dsigma0 = 7/16*Moment0.*(1./T_duration/k./Vs).^3/1e6; 
% stress drop estimated directly from STF is of great uncertainty
%Dsigma = median(Dsigma0); 
Dsigma = 1; % use the constant stressdrop results from Allmann and Shearer 2009

Shear_modulus=Vs.^2.*Rho;

radiation_efficiency_stressdrop=2*Shear_modulus.*total_RE_all_scaled./Dsigma0/1e6;
radiation_efficiency=2*Shear_modulus.*total_RE_all_scaled./Dsigma/1e6;  

%%
MinPeakProminence_TH=0.1;

% group colormap
GROUP_CMP=[32 100 154;
    60 173 162;
    245 212 93;
    238 86 60]/255;

% reorganize the group label (cutoff complexity only to 4 (>=4 -> 4) groups)
Complexity_number=group_label;
Complexity_number(Complexity_number>4) = 4;
Complexity_label = {'G1','G2','G3','G4'};

%% ===================== Reproducing Figure 1 Left half =================================
f1=figure(1);
%% first plot the distribution of clustered STF
subplot(1,3,[3])
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


%% Go through event and Plotting the stretched results
stretched_reference_STF = dtw_stretched_stfs.stretched_reference_STF;
stretched_STF = dtw_stretched_stfs.stretched_STF;
reference_STF_indice = dtw_stretched_stfs.reference_STF_indice;

for ic=1:length(Label_event)
    I_cluster=str2double(cell2mat(Label_event(ic))); % value of this cluster label
    % line color index
    event_indx_cluster=find(cluster_label==I_cluster); % event index within each cluster
    
    % get the center events ready (peak location)
    
    [STF_PKS,STF_PKS_LOCS]=findpeaks(stretched_reference_STF(I_cluster,:),...
        'MinPeakProminence',MinPeakProminence_TH*max(stretched_reference_STF(I_cluster,:)));
    
    GROUP_NUMBER = Complexity_number(reference_STF_indice(I_cluster));
    GROUP_LABEL = Complexity_label{GROUP_NUMBER};
    GROUP_color = GROUP_CMP(GROUP_NUMBER,:);
    
    %% Plotting stretched STFs
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
    
    %% Plotting unstretched STFs (normalized)
    Matrix_stf = normalized_series;
    subplot(1,4,1)
    MAX_temp2=max(max(Matrix_stf(event_indx_cluster,:)));
    plot((1:N_pts),Matrix_stf(event_indx_cluster,:)'/MAX_temp2+1-ic,'-','LineWidth',1,'Color',[0.5 0.5 0.5 0.7]);
    hold on
    %max(max(stretched_reference_STF(I_cluster,:)))
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
  [output_dir '/Fig1_left_SCARDEC_stretched_STFs.pdf'])
%keyboard
figure(1)
%clf

%% ===================== Reproducing Figure 2 Source parameter correlation =================================
% check the depth and focal mechanism distribution in each group
DEPTH_TH=80;

N_bins=20;
for icpx=1:4
    figure(65)
    % -------------- (a) Depth correlation -----------------
    subplot(5,2,[1.2 1.8])
    hold on
    
    I=find(Complexity_number==icpx & Depth <=DEPTH_TH);
    % Group depth distribution
    [Prob_density_depth,BinEdges_depth] = histcounts(Depth(I),linspace(0,80,N_bins));
    Prob_density_depth=Prob_density_depth/sum(Prob_density_depth);
    plot(BinEdges_depth(1:(end-1)),Prob_density_depth,'Color',GROUP_CMP(icpx,:),'LineWidth',4)
    xlabel('Depth (km)')
    ylabel('Probablity Freq.')
    title('Depth')
    grid on
    
end
legend('Group 1','Group 2','Group 3','Group 4','FontName','Times')

for icpx=1:4
    figure(65)
    % -------------- (b) Focal mechanism correlation -----------------
    subplot(5,2,[3.2 3.8])
    hold on
    
    I=find(Complexity_number==icpx & Depth <=DEPTH_TH);
    % Group focal mechanism distribution
    [Prob_density_FM,BinEdges_FM] = histcounts(fm_vector(I),linspace(-1,1,N_bins));
    Prob_density_FM=Prob_density_FM/sum(Prob_density_FM);
    plot(BinEdges_FM(1:(end-1)),Prob_density_FM,'Color',GROUP_CMP(icpx,:),'LineWidth',4)
    
    AAA=subplot(5,2,[3.2 3.8]);
    AAA.XTick=-1:0.5:1;
    AAA.XTickLabel={'-1 (Normal)','-0.5','0 (Strike Slip)','0.5','1 (Thrust)'};
    
    ylabel('Probablity Freq.')
    title('Focal mechanism')
    grid on
end

for icpx=1:4
    figure(65)
    % -------------- (c) Scaled energy correlation -----------------
    subplot(5,2,[5.2 5.8])
    hold on
    
    I=find(Complexity_number==icpx & Depth <=DEPTH_TH);
    % Group focal mechanism distribution
    [Prob_density_DS,BinEdges_DS] = histcounts(total_RE_all_scaled(I),logspace(-8,-4,N_bins));
    %[Prob_density_DS,BinEdges_DS] = histcounts(Dsigma(I),[0:10:200]);
    Prob_density_DS=Prob_density_DS/sum(Prob_density_DS);
    plot(BinEdges_DS(1:(end-1)),Prob_density_DS,'Color',GROUP_CMP(icpx,:),'LineWidth',4)
    set(gca,'XScale','log')
    xlabel('E_R/M_0','Interpreter','tex','FontName','Times')
    ylabel('Probablity Freq.')
    title('Scaled energy')
    grid on
end

% -------------- (d) Duration correlation -----------------
subplot(5,2,[7.2 9.8])
for II=1:4
    KK=find(Complexity_number==II);
    loglog(Moment0(KK),T_duration(KK),'ko','MarkerSize',8,'MarkerFaceColor',GROUP_CMP(II,:),'LineWidth',0.5)
    alpha(0.5)
    hold on
end
xlim([0 10*max(Moment0)])
ylim([0 2*max(T_duration)])
xlabel('Moment (Nm)')
ylabel('Duration (s)')
legend('Group1','Group2','Group3','Group4','Location','northwest','FontName','Times')
grid on
title('SCARDEC')

f65=figure(65);
f65.Position=[152 0 945 1344];
f65.PaperSize=f65.Position(3:4);
print('-dpdf','-painters',[output_dir '/Fig2_Source_parameter_distribution.pdf'])

%%find two events with close location, magnitude, but different complexity
f66=figure(66);
event_loc = [Lat Lon Depth Complexity_number]; % lat, lon, depth
I_lat = event_loc(:,1)>=-18 & event_loc(:,1)<=-16;
I_lon = event_loc(:,2)>=160 & event_loc(:,2)<=170;
I_depth = event_loc(:,3)<=30 & event_loc(:,3)>=20;
I_find = find(I_lat & I_lon & I_depth);
I_find = I_find([7 11]) % 7 and 11 are the chosen indice
event_find = event_loc(I_find,:)

T_duration(I_find)
T1=linspace(0,T_duration(I_find(1)),100);
T2=linspace(0,T_duration(I_find(2)),100);



subplot(211)
plot(T1,All_STFs(I_find(1),:),'color',GROUP_CMP(Complexity_number(I_find(1)),:));
hold on
plot(T2,All_STFs(I_find(2),:),'color',GROUP_CMP(Complexity_number(I_find(2)),:));
ylim([0 1e18])
xlabel('Time (s)');
ylabel('Moment rate (Nm/s)');

subplot(212)
bb(FocalMechanism(I_find(1),1:3),1,1,5,0,GROUP_CMP(Complexity_number(I_find(1)),:));
text(1,-5,'STF1','FontSize',20)
hold on
bb(FocalMechanism(I_find(2),1:3),20,1,5,0,GROUP_CMP(Complexity_number(I_find(2)),:))
text(20,-5,'STF2','FontSize',20)
axis equal
xxx= gca;
xxx.Visible ='off';
f66.Position = [68 50 656 1006];
f66.PaperSize = f66.Position(3:4);
print('-dpdf','-painters',[output_dir '/Fig2_subpanel_colocated_events.pdf'])



%% ===================== Reproducing Figure S4 =================================
% Comparison between number of peaks in the center event and each
% stretched individual event, also show that the magnitude correlation is
% not clear

f5=figure(5);
for icpx=1:4
    subplot(2,2,icpx)
    I=find(Complexity_number==icpx);
    N_peak_original=zeros(size(I));
    Mw_original=zeros(size(I));
    for i_event=1:length(I)
        [temp_PKS,temp_LOCs]=findpeaks(Matrix_stf(I(i_event),:),...
            'MinPeakProminence',MinPeakProminence_TH*max(Matrix_stf(I(i_event),:)));
        N_peak_original(i_event)=length(temp_PKS);
        Mw_original(i_event)=Mw(I(i_event));

    end
    Individual_PKS{icpx,1}=N_peak_original;
    Individual_PKS{icpx,2}=Mw_original;
    
    [n,c] = hist3([N_peak_original, Mw_original],{[0.5:1:6.5] [5.5:0.5:9.5]});
    [Complexity_grid,subevent_grid]=meshgrid(c{1},c{2});
    Frequency=n'/sum(n(:));
    Frequency_scatter=interp2(Complexity_grid,subevent_grid,Frequency,N_peak_original,Mw_original);
    
    plot([icpx icpx],[5.5 9.5],'-','LineWidth',30,'Color',[1 0 0 0.3])
    hold on
    scatter(N_peak_original,Mw_original,15*Mw(I),Frequency_scatter,'filled','LineWidth',0.5)
    caxis([0 0.2])
    colorbar
    xlabel('Number of prominent peaks in unstretched STFs')
    ylabel('Magnitude')
    xlim([0 7])
    ylim([5.5 9.5])
    grid on
    title(['Group ' num2str(icpx) '(' num2str(length(I)) ')'],'Color',GROUP_CMP(icpx,:)+[-0.1 -0.2 -0.2])
end

f5.Position=[50 50 1400 1000];
f5.PaperSize=f5.Position(3:4);
print('-dpdf',[output_dir '/FigS4_Peak_number_Magnitude.pdf'])


%% 2) Comparison between number of peaks in the center event and Gaussian subevent number
fh=figure(6);

for icpx=1:4
    
    I=find(Complexity_number==icpx);
    
    [n,c] = hist3([Complexity_number(I), Mw(I)],{[icpx+[-0.5 0]] [5.5:0.5:9.5]});
    [Complexity_grid,Mw_grid]=meshgrid(c{1},c{2});
    Frequency=n'/sum(n(:));
    
    Frequency_scatter=interp2(Complexity_grid,Mw_grid,Frequency,icpx,Mw(I));
    
    scatter(Complexity_number(I),Mw(I),500,Frequency_scatter,'s','filled')
    hold on
    box_X = [icpx-0.5 icpx+0.5 icpx+0.5 icpx-0.5 icpx-0.5];
    box_Y = [0 0 27 27 0];
    plot(box_X,box_Y,'-','LineWidth',3,'Color',GROUP_CMP(icpx,:))
    
end

xlim([0 5])
ylim([5.5 10])
colormap(jet)
caxis([0 0.3])
colorbar
AA=gca;
AA.XTick = 1:4;
AA.XTickLabel = {'Group 1','Group 2','Group 3','Group 4'};
AA.YTick = 6:1:9;
ylabel('Moment magnitude (Mw)')
print('-dpdf',[output_dir '/FigSX_Peak_number_vs_magnitude.pdf'])

%% ===================== Reproducing Figure S5 Stressdrop correlation =================================
fh=figure(8);
subplot(3,4,[1.5 3.5])
scatter(fm_vector,Dsigma0,50,[1 0 0],...
    'MarkerFaceColor','r','MarkerEdgeColor',[0.5 0.5 0.5],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
hold on
plot(fm_vector, Dsigma*ones(size(fm_vector)),'-b')


% possible check within each fptype bin...
d_win=0.1;
fm_vector_win1=[-1:d_win:(1-d_win)];
fm_vector_win2=fm_vector_win1+d_win;
fm_center=fm_vector_win1+d_win/2;
stress_drop_median=zeros(size(fm_center));
std_median=zeros(size(fm_center));

for I_win=1:length(fm_vector_win1);
    II_win=find(fm_vector>=fm_vector_win1(I_win) & fm_vector<fm_vector_win2(I_win));
    stress_drop_median(I_win)=median(Dsigma0(II_win));
    
    std_median(I_win)=std(log10(Dsigma0(II_win)));
    
    plot(fm_center(I_win)*[1 1],[10.^(-std_median(I_win)) 10.^std_median(I_win)],'-k','LineWidth',2);
end



plot(fm_center,stress_drop_median,'ks','MarkerFaceColor','c','MarkerSize',15)

set(gca,'YScale','log')
ylim([1e-3,1e3])
xlabel('fptype')
ylabel('Stress drop (MPa)')
legend('Individual event','Constant 1 MPa','Median stressdrop','FontName','Times','FontSize',15)
set(gca,'XTick',-1:0.5:1,'XTickLabel',{'-1 (Normal)','-0.5','0 (Strike Slip)','0.5','1 (Thrust)'})
grid on

for icpx=1:4 % stress drop
    subplot(3,4,[5 6])
    hold on
    
    I=find(Complexity_number==icpx);
    % Group focal mechanism distribution
    [Prob_density_Dsigma,BinEdges_Dsigma] = histcounts(Dsigma0(I),logspace(-2,2,N_bins));
    %[Prob_density_DS,BinEdges_DS] = histcounts(Dsigma(I),[0:10:200]);
    Prob_density_Dsigma=Prob_density_Dsigma/sum(Prob_density_Dsigma);
    plot(BinEdges_Dsigma(1:(end-1)),Prob_density_Dsigma,'Color',GROUP_CMP(icpx,:),'LineWidth',4)
    set(gca,'XScale','log')
    xlabel('\Delta\sigma_M (MPa)','Interpreter','tex','FontName','Times')
    ylabel('Probablity Freq.')
    title('Stress drop')
    grid on
end

% use individual stress drop
for icpx=1:4 % stress drop
    subplot(3,4,[9 10])
    hold on
    
    I=find(Complexity_number==icpx);
    % Group focal mechanism distribution
    [Prob_density_Dsigma,BinEdges_Dsigma] = histcounts(radiation_efficiency_stressdrop(I),logspace(-2,2,N_bins));
    %[Prob_density_DS,BinEdges_DS] = histcounts(Dsigma(I),[0:10:200]);
    Prob_density_Dsigma=Prob_density_Dsigma/sum(Prob_density_Dsigma);
    plot(BinEdges_Dsigma(1:(end-1)),Prob_density_Dsigma,'Color',GROUP_CMP(icpx,:),'LineWidth',4)
    set(gca,'XScale','log')
    xlabel('\eta_R=2\muE_R/M_0/\Delta\sigma_M','Interpreter','tex','FontName','Times')
    ylabel('Probablity Freq.')
    title('Radiation ratio from individual event')
    grid on
end

% assuming constant stress drop of 1MPa
radiation_efficiency_stressdrop=2*Shear_modulus.*total_RE_all_scaled./1e6;
for icpx=1:4 % stress drop
    subplot(3,4,[11 12])
    hold on
    
    I=find(Complexity_number==icpx);
    % Group focal mechanism distribution
    [Prob_density_Dsigma,BinEdges_Dsigma] = histcounts(radiation_efficiency_stressdrop(I),logspace(-2,2,N_bins));
    %[Prob_density_DS,BinEdges_DS] = histcounts(Dsigma(I),[0:10:200]);
    Prob_density_Dsigma=Prob_density_Dsigma/sum(Prob_density_Dsigma);
    plot(BinEdges_Dsigma(1:(end-1)),Prob_density_Dsigma,'Color',GROUP_CMP(icpx,:),'LineWidth',4)
    set(gca,'XScale','log')
    xlabel('\eta_R=2\muE_R/M_0/(1MPa)','Interpreter','tex','FontName','Times')
    ylabel('Probablity Freq.')
    title('Radiation ratio from constant stress drop')
    grid on
end

% Strain drop
strain_drop = Dsigma0./Shear_modulus;
for icpx=1:4 % stress drop
    subplot(3,4,[7 8])
    hold on
    
    I=find(Complexity_number==icpx);
    % Group focal mechanism distribution
    [Prob_density_Dsigma,BinEdges_Dsigma] = histcounts(strain_drop(I),logspace(-13,-9,N_bins));
    %[Prob_density_DS,BinEdges_DS] = histcounts(Dsigma(I),[0:10:200]);
    Prob_density_Dsigma=Prob_density_Dsigma/sum(Prob_density_Dsigma);
    plot(BinEdges_Dsigma(1:(end-1)),Prob_density_Dsigma,'Color',GROUP_CMP(icpx,:),'LineWidth',4)
    set(gca,'XScale','log')
    xlabel('\Delta\epsilon=\Delta\sigma_M /\mu','Interpreter','tex','FontName','Times')
    ylabel('Probablity Freq.')
    title('Strain drop')
    grid on
end


fh.Position=[50 50 1400 1000];
fh.PaperSize=fh.Position(3:4);
print('-dpdf',[output_dir '/FigS5_Stressdrop_vs_group.pdf'],'-painters')


%% Further look at the peak distributions within each group
%% 1) Comparison between number of peaks in the center event and each unstretched individual event
for icpx=1:4
    f4=figure(4);
    subplot(2,2,icpx)
    I=find(Complexity_number==icpx);
    N_peak_original=zeros(size(I));
    N_peak_stretch=zeros(size(I));
    for i_event=1:length(I)
        [temp_PKS,temp_LOCs]=findpeaks(Matrix_stf(I(i_event),:),...
            'MinPeakProminence',MinPeakProminence_TH*max(Matrix_stf(I(i_event),:)));
        N_peak_original(i_event)=length(temp_PKS);
        
        [temp_PKS,temp_LOCs]=findpeaks(stretched_STF(I(i_event),:),...
            'MinPeakProminence',MinPeakProminence_TH*max(stretched_STF(I(i_event),:)));
        N_peak_stretch(i_event)=length(temp_PKS);
    end
    Individual_PKS{icpx}=N_peak_original;
    Individual_stretched_PKS{icpx}=N_peak_stretch;
    
    
    histogram(N_peak_original,[0.5:1:6.5]);
    hold on
    histogram(N_peak_stretch,[0.5:1:6.5]);
    xlabel('Peaks number')
    lg=legend('Original','Stretched');
    lg.FontName='Times';
    title(['Group ' num2str(icpx) '(' num2str(length(I)) ')'],'Color',GROUP_CMP(icpx,:)+[-0.1 -0.2 -0.2])
    
end

f4.Position=[50 50 1400 1000]
f4.PaperSize=f4.Position(3:4);
print('-dpdf',[output_dir '/FigSX_Peak_number_vs_unstretched_stf.pdf'])

%% 2) Comparison between number of peaks in the center event and Gaussian subevent number
fh=figure(9);

for icpx=1:4
    
    I=find(Complexity_number==icpx);
    
    [n,c] = hist3([Complexity_number(I), gaussian_subevent_number(I)],{[icpx+[-0.5 0]] [0:1:25]});
    [Complexity_grid,subevent_grid]=meshgrid(c{1},c{2});
    Frequency=n'/sum(n(:));
    
    Frequency_scatter=interp2(Complexity_grid,subevent_grid,Frequency,icpx,gaussian_subevent_number(I));
    
    scatter(Complexity_number(I),gaussian_subevent_number(I),500,Frequency_scatter,'s','filled')
    hold on
    box_X = [icpx-0.5 icpx+0.5 icpx+0.5 icpx-0.5 icpx-0.5];
    box_Y = [0 0 27 27 0];
    plot(box_X,box_Y,'-','LineWidth',3,'Color',GROUP_CMP(icpx,:))
    
end

xlim([0 5])
ylim([0 28])
colormap(jet)
caxis([0 0.3])
colorbar
AA=gca;
AA.XTick = 1:4;
AA.XTickLabel = {'Group 1','Group 2','Group 3','Group 4'};
AA.YTick = 0:5:25;
ylabel('Number of Gaussian-subevents')
print('-dpdf',[output_dir '/FigSX_Peak_number_vs_subevent_number.pdf'])

%% Function 
function adding_text(N_pt,ic,GROUP_color,GROUP_LABEL,event_indx_cluster)
text(0.99*(N_pt),1-ic+0.5,['\color[rgb]{' num2str(GROUP_color) '}' GROUP_LABEL '(' num2str(length(event_indx_cluster)) ')'],...
    'FontSize',18,'HorizontalAlignment','left','Interpreter','Tex')
end





