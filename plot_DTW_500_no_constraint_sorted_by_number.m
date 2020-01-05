MinPeakProminence_TH=0.1;

% do either maxima scaling or area scaling
Scaling_Index_plot=1; % 0: stf x RE,  1: stf only
if Scaling_Index_plot==0
    txt_tail=['_DTW_area_scaling_median'];
elseif Scaling_Index_plot==1
    txt_tail=['_DTW_500_area_scaling_median_STF_only'];
    
    stretched_reference_RE=zeros(20,100)*nan;
end

CMP=[32 100 154;
    60 173 162;
    245 212 93;
    238 86 60]/255;
%% load STFs
if Scaling_Index_plot==0
    load(['stfre_info.mat']);
else
    load([FILE_name '/clustered_STFs_x_SAF_PW_fix_phase' txt_tail '_DTW_info.mat'])
    load(['All_stfs_500.mat']);
    fm_vector=focal_mechanism_vector(All_headers(:,[14 17])); % focal mechanism of each events
end

Complexity=cell(size(DTW_labels)); % to store the complexity labels
Num_complexity = zeros(size(DTW_labels)); % to store the complexity labels (G1:0; G2:1; G3:2; G4:3)
N_pt=size(stretched_STF,2);

f1=figure(1);
subplot(1,3,[3])
H=histogram(DTW_labels,'DisplayOrder','ascend');
H.Orientation='horizontal';
H.Normalization='probability';
H.FaceColor='b';
set(gca, 'box','off','YTickLabel',[],'YTick',[],'XScale','linear')
xlim([0 0.4])
grid on
xlabel('Number of event')
title('Distributions')

[Num_of_event,Label_event]=histcounts(DTW_labels);

% Resort the label of event
[Num_of_event,Indx_sort]=sort(Num_of_event,'descend');
Label_event=Label_event(Indx_sort);


%% Go through event and Plotting the stretched results
for ic=1:length(Label_event)
    I_cluster=str2double(cell2mat(Label_event(ic))); % value of this cluster label
    % line color index
    
    event_indx_cluster=find(DTW_labels==Label_event(ic)); % event index within each cluster
    
    %-------------------% For temporal use
    if Scaling_Index_plot==1
        stretched_reference_RE(I_cluster,:)=mean(stretched_RE(event_indx_cluster,:),1);
    end
    %-------------------%
    
    [STF_PKS,STF_PKS_LOCS]=findpeaks(stretched_reference_STF(I_cluster,:),...
        'MinPeakProminence',MinPeakProminence_TH*max(stretched_reference_STF(I_cluster,:)));
    [GROUP_LABEL,GROUP_NUMBER,GROUP_color]=group_and_color([1,2,3],length(STF_PKS),CMP);
    
    if Scaling_Index_plot~=1
        [RE_PKS,RE_PKS_LOCS]=findpeaks(stretched_reference_RE(I_cluster,:),...
            'MinPeakProminence',MinPeakProminence_TH*max(stretched_reference_RE(I_cluster,:)));
        [GROUP_LABEL,GROUP_NUMBER,GROUP_color]=group_and_color([3,12,27],length(STF_PKS)*length(RE_PKS),CMP);
    end
    
    % complexity labels
    Complexity(event_indx_cluster)={GROUP_LABEL};
    Num_complexity(event_indx_cluster)=GROUP_NUMBER;
    
    %Plotting distributions and STFs
    MAX_temp1=max(max(stretched_STF(event_indx_cluster,:)))/1;
    subplot(1,12,[5 7])
    plot(1:N_pt,stretched_STF(event_indx_cluster,:)'/MAX_temp1+1-ic,'-','LineWidth',1,'Color',[0.5 0.5 0.5 0.7]);
    hold on
    plot(1:N_pt,stretched_reference_STF(I_cluster,:)/MAX_temp1+1-ic,'-k','LineWidth',2.5);
    plot(STF_PKS_LOCS,STF_PKS/MAX_temp1+1-ic,'o','MarkerSize',8,'Color',GROUP_color,'MarkerFaceColor',GROUP_color);
    
    title('Stretched STFs')
    set(gca, 'box','off','XColor','w','YColor','w','TickDir','out')
    adding_text(N_pt,ic,GROUP_color,GROUP_LABEL,event_indx_cluster)
    ylim([-19 1])
    
    
    subplot(1,4,1)
    if Scaling_Index_plot==1 % STF cluster only, RE is stretched based on STF matching
        %MAX_temp2=max(max(stretched_RE(event_indx_cluster,:)))/1;
        %plot((1:N_pt),stretched_RE(event_indx_cluster,:)'/MAX_temp2+1-ic,'-','LineWidth',1,'Color',[0.5 0.5 0.5 0.7]);
        
        MAX_temp2=max(max(Matrix_stf(event_indx_cluster,:)))/1;
        plot((1:N_pt),Matrix_stf(event_indx_cluster,:)'/MAX_temp2+1-ic,'-','LineWidth',1,'Color',[0.5 0.5 0.5 0.7]);
        hold on
        plot((1:N_pt),stretched_reference_STF(I_cluster,:)/MAX_temp2+1-ic,'-k','LineWidth',2.5);
        adding_text(N_pt,ic,GROUP_color,GROUP_LABEL,event_indx_cluster)
        titleTEXT='Unstretched STFs';
    else % STF x RE cluster
        MAX_temp2=max(max(stretched_RE(event_indx_cluster,:)))/1.1;
        plot((1:N_pt),stretched_RE(event_indx_cluster,:)'/MAX_temp2+1-ic,'-','LineWidth',1,'Color',[0.5 0.5 0.5 0.7]);
        hold on
        plot((1:N_pt),stretched_reference_RE(I_cluster,:)/MAX_temp2+1-ic,'-k','LineWidth',2.5);
        plot(RE_PKS_LOCS,RE_PKS/MAX_temp2+1-ic,'o','MarkerSize',8,'Color',GROUP_color,'MarkerFaceColor',GROUP_color);
        adding_text(N_pt,ic,GROUP_color,GROUP_LABEL,event_indx_cluster)
        titleTEXT='Stretched REs';
    end
    
    ylim([-19 1])
    title(titleTEXT);
    set(gca, 'box','off','XColor','w','YColor','w','TickDir','out')
    
end
f1.Position=[50 50 800 927];
f1.PaperSize=f1.Position(3:4);

%keyboard
print('-dpdf','-painters',...
   ['SCARDEC' txt_tail '_prominence_' num2str(MinPeakProminence_TH) '.pdf'])
%keyboard
figure(1)
%clf


%% Further grouping
Complexity_number=zeros(size(Complexity));
Complexity_number(categorical(Complexity)=='G1')=1;
Complexity_number(categorical(Complexity)=='G2')=2;
Complexity_number(categorical(Complexity)=='G3')=3;
Complexity_number(categorical(Complexity)=='G4')=4;

DEPTH_TH=80;
f2=figure(2);
f2.Position=[59 81 1061 852];
f2.PaperSize=f2.Position(3:4);
for icpx=1:4
    figure(2)
    subplot(4,2,2*(icpx-1)+1)
    yyaxis left
    I=find(Complexity_number==icpx & Depth <=DEPTH_TH);
    scatter(Depth(I),Mg(I),15*Mg(I),Complexity_number(I),'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor',CMP(icpx,:),'LineWidth',0.5)
    set(gca,'YColor','k')
    ylabel('Mw')
    xlim([0 DEPTH_TH])
    ylim([5.5 9.5])
    alpha(0.8)
    
    % subplot(4,2,2*icpx)
    % h=histogram(Depth(I),[0:5:110]);
    % h.FaceColor=CMP(icpx,:);
    % h.FaceAlpha=1;
    % h.Normalization='probability';
    % xlim([0 100])
    % ylim([0 0.5])
    
    
    [Prob_density,BinEdges] = histcounts(Depth(I),[0:5:110]);
    Prob_density=Prob_density/sum(Prob_density);
    
    yyaxis right
    plot(BinEdges(1:(end-1)),Prob_density,'Color','m')
    ylim([0 0.5])
    set(gca,'YColor','m')
    ylabel('Frequency')
    
    title(['Group ' num2str(icpx) '(' num2str(length(I)) ')'],'Color',CMP(icpx,:)+[-0.1 -0.2 -0.2])
    
    f3=figure(3);
    subplot(4,2,2*(icpx-1)+1)
    yyaxis left
    [n,c] = hist3([Depth(I), Mg(I)],{[0:5:110] [5.5:0.5:9.5]});
    [DEP_grid,MG_grid]=meshgrid(c{1},c{2});
    Frequency=n'/sum(n(:));
    Frequency_scatter=interp2(DEP_grid,MG_grid,Frequency,Depth(I),Mg(I));
    
    scatter(Depth(I),Mg(I),15*Mg(I),Frequency_scatter,'filled','LineWidth',0.5)
    alpha(0.8)
    colormap(jet)
    %colorbar
    caxis([0 0.3])
    ylabel('Mw')
    xlim([0 100]);
    ylim([5.5 9.5])
    alpha(0.8)
    title(['Group ' num2str(icpx) '(' num2str(length(I)) ')'],'Color',CMP(icpx,:)+[-0.1 -0.2 -0.2])
    
    [Prob_density,BinEdges] = histcounts(Depth(I),[0:5:110]);
    Prob_density=Prob_density/sum(Prob_density);
    
    yyaxis right
    plot(BinEdges(1:(end-1)),Prob_density,'Color','m')
    ylim([0 0.5])
    set(gca,'YColor','m')
    ylabel('Frequency')
end

f3.Position=[59 81 1061 852];
f3.PaperSize=f3.Position(3:4);
figure(3)
%print('-dpdf',['SARDEC_Frequency_' txt_tail '_prominence_' num2str(MinPeakProminence_TH) '.pdf'],'-painters')

figure(2)
subplot(422)
yyaxis left
I=find(Depth > DEPTH_TH);
scatter(Depth(I),Mg(I),15*Mg(I),Complexity_number(I),'filled','MarkerEdgeColor','k','LineWidth',0.5)
xlim([DEPTH_TH 700])
ylim([5.5 9.5])
grid on
colormap(CMP)
alpha(0.8)
set(gca,'YColor','k')

[Prob_density,BinEdges] = histcounts(Depth(I),[0:5:110]);
Prob_density=Prob_density/sum(Prob_density);

yyaxis right
plot(BinEdges(1:(end-1)),Prob_density,'Color','m')
ylim([0 0.5])
set(gca,'YColor','m')
ylabel('Density')

title(['Deep events' '(' num2str(length(I)) ')'])

subplot(4,2,[4 8])
HH=histogram(categorical(Complexity));
HH.Normalization='probability';
title('Distribution')
ylabel('Density')
ylim([0 1])

%print('-dpdf',['SARDEC_depth_' txt_tail '_prominence_' num2str(MinPeakProminence_TH) '.pdf'],'-painters')


%% check the depth and focal mechanism distribution in each group
N_bins=20;
for icpx=1:4
    figure(65)
    subplot(221)
    hold on
    
    I=find(Complexity_number==icpx & Depth <=DEPTH_TH);
    % Group depth distribution
    [Prob_density_depth,BinEdges_depth] = histcounts(Depth(I),linspace(0,80,N_bins));
    Prob_density_depth=Prob_density_depth/sum(Prob_density_depth);
    plot(BinEdges_depth(1:(end-1)),Prob_density_depth,'Color',CMP(icpx,:),'LineWidth',4)
    xlabel('Depth (km)')
    ylabel('Probablity Freq.')
    title('Depth')
    grid on
    
end
legend('Group 1','Group 2','Group 3','Group 4','FontName','Times')

for icpx=1:4
    figure(65)
    subplot(222)
    hold on
    
    I=find(Complexity_number==icpx & Depth <=DEPTH_TH);
    % Group focal mechanism distribution
    [Prob_density_FM,BinEdges_FM] = histcounts(fm_vector(I),linspace(-1,1,N_bins));
    Prob_density_FM=Prob_density_FM/sum(Prob_density_FM);
    plot(BinEdges_FM(1:(end-1)),Prob_density_FM,'Color',CMP(icpx,:),'LineWidth',4)
    
    AAA=subplot(222);
    AAA.XTick=-1:0.5:1;
    AAA.XTickLabel={'-1 (Normal)','-0.5','0 (Strike Slip)','0.5','1 (Thrust)'};
    
    ylabel('Probablity Freq.')
    title('Focal mechanism')
    grid on
end

for icpx=1:4
    figure(65)
    subplot(223)
    hold on
    
    I=find(Complexity_number==icpx & Depth <=DEPTH_TH);
    % Group focal mechanism distribution
    [Prob_density_DS,BinEdges_DS] = histcounts(total_RE_all_scaled(I),logspace(-8,-4,N_bins));
    %[Prob_density_DS,BinEdges_DS] = histcounts(Dsigma(I),[0:10:200]);
    Prob_density_DS=Prob_density_DS/sum(Prob_density_DS);
    plot(BinEdges_DS(1:(end-1)),Prob_density_DS,'Color',CMP(icpx,:),'LineWidth',4)
    set(gca,'XScale','log')
    xlabel('E_R/M_0','Interpreter','tex','FontName','Times')
    ylabel('Probablity Freq.')
    title('Scaled energy')
    grid on
end

for icpx=1:4 % radiation efficiency
    figure(65)
    subplot(224)
    hold on
    
    I=find(Complexity_number==icpx);
    % Group focal mechanism distribution
    [Prob_density_eta,BinEdges_eta] = histcounts(radiation_efficiency(I),logspace(-3,1,N_bins));
    %[Prob_density_DS,BinEdges_DS] = histcounts(Dsigma(I),[0:10:200]);
    Prob_density_eta=Prob_density_eta/sum(Prob_density_eta);
    plot(BinEdges_eta(1:(end-1)),Prob_density_eta,'Color',CMP(icpx,:),'LineWidth',4)
    set(gca,'XScale','log')
    xlabel('\muE_R/M_0/\Delta\sigma','Interpreter','tex','FontName','Times')
    ylabel('Probablity Freq.')
    title('Radiation efficiency')
    grid on
end


% for icpx=1:4
%     figure(65)
%     subplot(224)
%     hold on
%     
%     I=find(Complexity_number==icpx & Depth <=DEPTH_TH);
%     % Group focal mechanism distribution
%     [Prob_density_RE,BinEdges_RE] = histcounts(log10(total_RE_all_scaled(I)),[-8:0.4:-4]);
%     %[Prob_density_DS,BinEdges_DS] = histcounts(Dsigma(I),[0:10:200]);
%     Prob_density_RE=Prob_density_RE/sum(Prob_density_RE);
%     plot(BinEdges_RE(1:(end-1)),Prob_density_RE,'Color',CMP(icpx,:),'LineWidth',4)
%     
%     xlabel('log_{10}(E_R/M0)','Interpreter','tex')
%     ylabel('Frequency of occurrence')
%     title('Scaled radiated energy')
% end

% subplot(221)
% legend('G1','G2','G3','G4')
% subplot(222)
% legend('G1','G2','G3','G4')
% subplot(223)
% legend('G1','G2','G3','G4')
% subplot(224)
% legend('G1','G2','G3','G4')

f65=figure(65);
f65.Position=[99 123 1524 645];
f65.PaperSize=f65.Position(3:4);
print('-dpdf','-painters','Source_parameter_distribution.pdf')

%% Further look at the peak distributions within each group
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
    title(['Group ' num2str(icpx) '(' num2str(length(I)) ')'],'Color',CMP(icpx,:)+[-0.1 -0.2 -0.2])
    
end

f4.Position=[50 50 1400 1000]
f4.PaperSize=f4.Position(3:4);


print('-dpdf',['Peak_number' txt_tail '.pdf'])

%%
f5=figure(5);
for icpx=1:4
    subplot(2,2,icpx)
    I=find(Complexity_number==icpx);
    N_peak_original=zeros(size(I));
    Mg_original=zeros(size(I));
    for i_event=1:length(I)
        [temp_PKS,temp_LOCs]=findpeaks(Matrix_stf(I(i_event),:),...
            'MinPeakProminence',MinPeakProminence_TH*max(Matrix_stf(I(i_event),:)));
        N_peak_original(i_event)=length(temp_PKS);
        Mg_original(i_event)=Mg(I(i_event));

    end
    Individual_PKS{icpx,1}=N_peak_original;
    Individual_PKS{icpx,2}=Mg_original;
    
    [n,c] = hist3([N_peak_original, Mg_original],{[0.5:1:6.5] [5.5:0.5:9.5]});
    [Peak_grid,Mg_grid]=meshgrid(c{1},c{2});
    Frequency=n'/sum(n(:));
    Frequency_scatter=interp2(Peak_grid,Mg_grid,Frequency,N_peak_original,Mg_original);
    
    plot([icpx icpx],[5.5 9.5],'-','LineWidth',30,'Color',[1 0 0 0.3])
    hold on
    scatter(N_peak_original,Mg_original,15*Mg(I),Frequency_scatter,'filled','LineWidth',0.5)
    caxis([0 0.2])
    colorbar
    xlabel('Peaks number')
    ylabel('Magnitude')
    xlim([0 7])
    ylim([5.5 9.5])
    grid on
    title(['Group ' num2str(icpx) '(' num2str(length(I)) ')'],'Color',CMP(icpx,:)+[-0.1 -0.2 -0.2])
end

f5.Position=[50 50 1400 1000]
f5.PaperSize=f5.Position(3:4);
print('-dpdf',['Peak_number_Magnitude' txt_tail '.pdf'])



%% Functions
function [GROUP_LABEL,GROUP_NUMBER,GROUP_color]=group_and_color(A,X,CMP)
% get the group complexity label and label color for each cluster

if X<=A(1)
    GROUP_color=CMP(1,:);
    GROUP_LABEL='G1';
    GROUP_NUMBER=0;
elseif X<=A(2)
    GROUP_color=CMP(2,:);
    GROUP_LABEL='G2';
    GROUP_NUMBER=1;
elseif X<=A(3)
    GROUP_color=CMP(3,:);
    GROUP_LABEL='G3';
    GROUP_NUMBER=2;
else
    GROUP_color=CMP(4,:);
    GROUP_LABEL='G4';
    GROUP_NUMBER=3;
end

end

function adding_text(N_pt,ic,GROUP_color,GROUP_LABEL,event_indx_cluster)
text(0.99*(N_pt),1-ic+0.5,['\color[rgb]{' num2str(GROUP_color) '}' GROUP_LABEL '(' num2str(length(event_indx_cluster)) ')'],...
    'FontSize',18,'HorizontalAlignment','left','Interpreter','Tex')
end




