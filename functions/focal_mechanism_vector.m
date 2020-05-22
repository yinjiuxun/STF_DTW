function fm_vector=focal_mechanism_vector(rake)
% 
% fm_vector=focal_mechanism_vector(rake)
% rake is Nx2 rake matrix in degree (rake1, rake2)
% vectorized focal mechanism mentioned by Shearer et al, 2006
% fm_vector parameter varies from -1 (normal) to 0 (strike slip) to 1 (thrust)
% and has the advantage of providing a single scalar value for characterizing the faulting type

rake1=rake(:,1);
rake2=rake(:,2);

II1=find(abs(rake1)>90);
rake1(II1)=(180-abs(rake1(II1))).*(rake1(II1)./abs(rake1(II1)));

II2=find(abs(rake2)>90);
rake2(II2)=(180-abs(rake2(II2))).*(rake2(II2)./abs(rake2(II2)));

clear II1 II2

fm_vector=zeros(size(rake,1),1);

II1=abs(rake1)<abs(rake2);
fm_vector(II1)=rake1(II1)/90;
II2=abs(rake1)>=abs(rake2);
fm_vector(II2)=rake2(II2)/90;


end


