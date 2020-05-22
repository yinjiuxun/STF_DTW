function Simulated_Data = read_simulation(STF_file,N_pts)
% Load the simulated STFs from stochastic simulations
% Simulated_Data = read_scardec_stf(STF_file,N_pts)
% 
% Simulated_Data is a strucutre:
%        All_STFs: STF matrix
%         All_REs: Radiated energy matrix
%     T_duration: Duration vector
%      Event_info: Information of the events

%%
temp_load = load(STF_file);

N_stf = length(temp_load.subevent);
% initialize the output variables
Simulated_Data = struct;
T_duration = zeros(N_stf,1);
All_STFs=zeros(N_stf,N_pts);
All_REs=zeros(N_stf,N_pts);
Event_info = struct;

% Event inforation
Event_info.Mg=log10([temp_load.subevent.MainMag_smooth]); % log10 "magnitude"
Event_info.X_extend=reshape([temp_load.subevent.X_extend],2,N_stf); % rupture extension in space
Event_info.Length=Event_info.X_extend(2,:)-Event_info.X_extend(1,:); % rupture length
Event_info.T_extend=reshape([temp_load.subevent.T_extend],2,N_stf); % time extension
Event_info.Duration=Event_info.T_extend(2,:)-Event_info.T_extend(1,:); % rupture duration
Event_info.moment_sim=[temp_load.subevent.moment_sim];

% energy calculated directly from energy partitioning during simulations
Event_info.total_energy=[temp_load.subevent.total_energy];
Event_info.fracture_energy=[temp_load.subevent.fracture_energy];
Event_info.radiated_energy=[temp_load.subevent.radiated_energy];

% different types of stressdrop
Event_info.Dsigma_average=[temp_load.subevent.Dsigma_average];
Event_info.Dsigma_energy=[temp_load.subevent.Dsigma_energy];
Event_info.Dsigma_duration=[temp_load.subevent.Dsigma_duration]; 
% for 2D simulation, stress drop estimated from duration is not correct due
% to the missing dimension of moment!



disp(['   ']);
disp(['Loading simualted STFs from: ' STF_file]);
disp(['   ']);
for istf=1:N_stf
    % Time and deltaT of each simulated STFs
    temp_T = temp_load.subevent(istf).T;
    dt=temp_T(2)-temp_T(1);
    
    % duration of events
    T_duration(istf) = Event_info.Duration(istf);
    
    % STFs and radiated energy
    temp_STF=temp_load.subevent(istf).stf;
    temp_RE=abs(diff(temp_STF)/dt).^2;
    
    % to scale and down sample both the STF and SAF
    temp_STF=temp_STF(temp_T>=(temp_load.subevent(istf).T_extend(1)-1) ...
        & temp_T<=(temp_load.subevent(istf).T_extend(2)+1));
    temp_RE=temp_RE(temp_T>=(temp_load.subevent(istf).T_extend(1)-1) ...
        & temp_T<=(temp_load.subevent(istf).T_extend(2))+1);
    
    temp_STF=interp1(1:length(temp_STF),temp_STF,linspace(1,length(temp_STF),N_pts));
    temp_STF(isnan(temp_STF))=0;
    %temp_STF_normal=temp_STF_normal/max(abs(temp_STF_normal));
    
    temp_RE=interp1(1:length(temp_RE),temp_RE,linspace(1,length(temp_RE),N_pts));
    temp_RE(isnan(temp_RE))=0;
     
    
    All_STFs(istf,:) = temp_STF;
    All_REs(istf,:) = temp_RE;

     
       
end

%% Assemble output structure
Simulated_Data.All_STFs = All_STFs;
Simulated_Data.All_REs = All_REs;
Simulated_Data.T_duration = T_duration;
Simulated_Data.Event_info = Event_info;






