function Scardec_Data = read_scardec_stf(STF_dir,N_pts)
% Load the SCARDEC STFs in the original form directly downloaded from website
% Scardec_Data = read_scardec_stf(STF_dir,N_pts)
% Scardec_Data is a strucutre:
%        All_STFs: STF matrix
%         All_REs: Radiated energy matrix
%     T_duration: Duration vector
%      Event_info: Information of the events

%%
FCTs=dir([STF_dir '/FCTs_*']);

% initialize the output variables
Scardec_Data = struct;
T_duration = zeros(length(FCTs),1);
All_STFs=zeros(length(FCTs),N_pts);
All_REs=zeros(length(FCTs),N_pts);
Event_info = struct;

% temp variables
All_headers_temp=zeros(length(FCTs),17);
gaussian_subevent_number = zeros(length(FCTs),1);


disp(['   ']);
disp(['Loading SCARDEC STFs from: ' STF_dir]);
disp(['   ']);
for ifcts=1:length(FCTs)
    
    if mod(ifcts,100)==0
        disp([int2str(ifcts) '/' int2str(length(FCTs)) ' loaded...']);
    end
    
    stf_file0=dir([FCTs(ifcts).folder '/' FCTs(ifcts).name '/fctmoysource_*']);
    stf_file=[stf_file0.folder '/' stf_file0.name];
    
    headers_temp=load_scardec_headers(stf_file);
    
    [T0, STF0] = load_scardec(stf_file);
      
    % radiated energy function
    RE=(diff(STF0)./diff(T0)).^2;
    T_RE=T0(1:(end-1))+diff(T0)/2;
    RE0=interp1(T_RE,RE,T0);
    RE0(isnan(RE0))=0;
    
    
    % determine the trucation time
    percent = 0.999;
    T_trunc = truncation_time (T0,STF0,percent);
    T_duration(ifcts)=T_trunc;
    
    
    % resample STF
    T_STF_resample = linspace(T0(1),T_trunc,N_pts);
    STF_temp=interp1(T0(T0<=T_trunc),STF0(T0<=T_trunc),T_STF_resample);
    RE_temp=interp1(T0(T0<=T_trunc),RE0(T0<=T_trunc),T_STF_resample);

    % remove possible nan
    STF_temp(isnan(STF_temp))=0;
    RE_temp(isnan(RE_temp))=0;
    
    All_STFs(ifcts,:)=STF_temp;
    All_REs(ifcts,:)=RE_temp;
    All_headers_temp(ifcts,:)=headers_temp;
    
    % decompose the STF into Gaussian subevents (Danre et al., 2019)
    subevent_threshold = 0.1;
    %disp(['Decomposing into Gaussian-shape subevents with threshold ' num2str(subevent_threshold) ' of maximum...'])
    GAUSSlet=fit_gaussian_peak(T_STF_resample,STF_temp,subevent_threshold);
    gaussian_subevent_number(ifcts) = length(GAUSSlet.amp);
       
end

%% STF source parameters from SCARDEC website
Event_info.Event_time = All_headers_temp(:,1:6);
Event_info.Lat = All_headers_temp(:,7);
Event_info.Lon = All_headers_temp(:,8);
Event_info.Depth = All_headers_temp(:,9);
Event_info.Mw = All_headers_temp(:,11);
Event_info.Moment = All_headers_temp(:,10);
Event_info.FocalMechanism = All_headers_temp(:,12:17);
Event_info.gaussian_subevent_number = gaussian_subevent_number;

%% Assemble output structure
Scardec_Data.All_STFs = All_STFs;
Scardec_Data.All_REs = All_REs;
Scardec_Data.T_duration = T_duration;
Scardec_Data.Event_info = Event_info;






