function GAUSSlet=fit_gaussian_peak(T,STF,threshold)



dt=T(2)-T(1);

RESIDUAL=STF;

STF_temp=STF;

%[PKS,LOCS]= findpeaks(STF_temp);%,'MinPeakProminence',0.1*max(STF));%,'MinPeakDistance',round(0.5/(dt)));

GAUSSlet.bestfit=zeros(size(T));

I_pks_temp=[1];
II0=0;
while ~isempty(I_pks_temp)
    
    [PKS,LOCS]= findpeaks(STF_temp);
    I_pks_temp=find(PKS>threshold*max(STF));
    
    if isempty(I_pks_temp)
        break
    else
        I_pks=I_pks_temp(1);
    end
    

    AMP0=STF_temp(LOCS(I_pks));
    MEAN0=T(LOCS(I_pks));
    SIGMA_test=[0.5:0.5:100];
    
    if I_pks==1
        Inx_fit=1:LOCS(I_pks);
    else
        %Inx_fit=round(LOCS(I_pks)-(LOCS(I_pks)-LOCS(I_pks-1))/2):LOCS(I_pks);
        Inx_fit=LOCS(I_pks-1):LOCS(I_pks);
    end
    
    DATA_STF=STF_temp(Inx_fit);
    DATA_T=T(Inx_fit);
    
    misfit=zeros(size(SIGMA_test)); % fitting misfit
    for ii=1:length(SIGMA_test)
        GAUSS_test=gauss_shape(DATA_T,MEAN0,SIGMA_test(ii),AMP0);
        
        misfit(ii)=norm(GAUSS_test-DATA_STF,2);
    end
    
    [~,I_best]=min(misfit);
    
    SIGMA_best=SIGMA_test(I_best);
    
    GAUSS_best=gauss_shape(T,MEAN0,SIGMA_best,AMP0);
    
    II0=II0+1;
    
%     figure(99)
%     subplot(2,1,1)
%     hold off
%     plot(T,STF_temp)
%     hold on
%     plot(T(LOCS),STF_temp(LOCS),'rx')
%     plot(T,GAUSS_best)
%     plot(T(Inx_fit),STF_temp(Inx_fit),':k')
        
    
    RESIDUAL=STF_temp-GAUSS_best;
    STF_temp=RESIDUAL;
    

%     subplot(2,1,2)
%     hold off
%     plot(T,RESIDUAL)
%     pause
    
    GAUSSlet.mean(II0)=MEAN0;
    GAUSSlet.sigma(II0)=SIGMA_best;
    GAUSSlet.amp(II0)=AMP0;    
    GAUSSlet.bestfit=GAUSSlet.bestfit+gauss_shape(T,MEAN0,SIGMA_best,AMP0);
    
end


function G=gauss_shape(x,mean,sigma,amp)

G=normpdf(x,mean,sigma);
G=G/max(G)*amp;






