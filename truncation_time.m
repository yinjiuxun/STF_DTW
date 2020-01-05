function T_trunc = truncation_time (T0,STF0,percent)
 % get the trucation time based on certain percent of total moment and
 % radiated energy
 
 % T_trunc = truncation_time (T0,STF0,percent)
 
 
 % getting RE 
 RE=(diff(STF0)./diff(T0)).^2;
 T_RE=T0(1:(end-1))+diff(T0)/2;
 
 RE0=interp1(T_RE,RE,T0);
 RE0(isnan(RE0))=0;
 
 total_RE=cumtrapz(T0,RE0);
 
 I_T_RE=find(total_RE<=percent*total_RE(end));
 T_RE_trunc=T0(I_T_RE(end)); % truncation time of radiated energy
 
 
% getting STF
 total_STF=cumtrapz(T0,STF0);
 
 I_T_STF=find(total_STF<=percent*total_STF(end));
 T_STF_trunc=T0(I_T_STF(end)); % truncation time of radiated energy

 
 T_trunc = (T_RE_trunc+T_STF_trunc)/2;

end