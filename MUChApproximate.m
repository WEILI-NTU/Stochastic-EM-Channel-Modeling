function [TheoChannel_Ac,TheoChannel_FF,TheoChannel_FF_ap]=MUChApproximate(par,temp_UETX_dis,termConstant,flagLoS,flagNLoS,users)
% This function computes the approximation of the channel capacity of MISO
%  for multi-users

for user_index=1:users.num
UETX_dis=par.UETX.distance(user_index,:);
constant=par.wavenumber*UETX_dis;
temp_direc_x= -temp_UETX_dis.X(user_index,:)./UETX_dis;
direc.XX=temp_direc_x.*temp_direc_x;

%% Accurate channel statistics (theo)

Exp_MISO_LoS_Ac=cos(constant)./(4*pi*UETX_dis)-sin(constant)./(4*pi*constant.*UETX_dis)...
    -cos(constant)./(4*pi*constant.^2.*UETX_dis)...
    +(3*cos(constant)./(4*pi*constant.^2.*UETX_dis)...
    +3*sin(constant)./(4*pi*constant.*UETX_dis)...
    -cos(constant)./(4*pi*UETX_dis)).*direc.XX;
Exp_MISO_NLoS_Ac=pi*par.Qfactor/(2*par.cavity*par.wavenumber^2)*(sin(constant)./constant...
    +cos(constant)./(constant.^2)-sin(constant)./(constant.^3));
Exp_MISO_Ac= flagLoS*termConstant*Exp_MISO_LoS_Ac+flagNLoS*termConstant*Exp_MISO_NLoS_Ac;
 
  
Var_MISO_Ac=9*pi^2*par.Qfactor^2/(128*par.cavity^2*par.wavenumber^4)*(4/3+...
    2*sin(2*constant)./(2*constant)+2*cos(2*constant)./(2*constant).^2  ...
    -2*sin(2*constant)./(2*constant).^3)-pi^2*par.Qfactor^2/(4*par.cavity^2 ...
    *par.wavenumber^4)*(sin(constant)./constant+constant.*cos(constant)./(constant.^3) ...
    -sin(constant)./constant.^3).^2;
Var_MISO_Ac=flagNLoS*Var_MISO_Ac;


TheoChannel_Ac(user_index,:)=Exp_MISO_Ac.*conj(Exp_MISO_Ac)+Var_MISO_Ac;


%% Approximate channel statistics in far-field
Exp_MISO_LoS_FF=cos(constant)./(4*pi*UETX_dis)*par.UETX.distance(1)^2./UETX_dis.^2;
Exp_MISO_NLoS_FF=pi*par.Qfactor/(2*par.cavity*par.wavenumber^2)*(sin(constant)./constant);
Exp_MISO_FF= flagLoS*termConstant*Exp_MISO_LoS_FF+flagNLoS*termConstant*Exp_MISO_NLoS_FF;

Var_MISO_FF=9*pi^2*par.Qfactor^2/(128*par.cavity^2*par.wavenumber^4)*(4/3+...
    2*sin(2*constant)./(2*constant))-pi^2*par.Qfactor^2/(4*par.cavity^2 ...
    *par.wavenumber^4)*(sin(constant)./constant).^2;
Var_MISO_FF=flagNLoS*Var_MISO_FF;

TheoChannel_FF(user_index,:)=Exp_MISO_FF.*conj(Exp_MISO_FF)+Var_MISO_FF; 


% sum(Exp_MISO_Ac)  
% sum(Exp_MISO_FF)
% sum(TheoChannel_Ac)
% sum(TheoChannel_FF)

%% Approximate channel statistics in far-field with negligible aperture
Exp_MISO_LoS_FF_ap=cos(constant(1))^2/(4*pi*UETX_dis(1))^2;
Exp_MISO_NLoS_FF_ap= par.Qfactor/(4*par.cavity*par.wavenumber) ...
    *(cos(constant(1))/constant(1))*(sin(constant(1))/constant(1)) ...
    +3*pi^2*par.Qfactor^2/(32*par.cavity^2*par.wavenumber^4) ...
    + 9*pi^2*par.Qfactor^2/(64*par.cavity^2*par.wavenumber^4) ...
    *(sin(2*constant(1))/(2*constant(1)));
Exp_MISO_FF_ap= flagLoS*termConstant*conj(termConstant) ...
    *Exp_MISO_LoS_FF_ap+flagNLoS*termConstant*conj(termConstant) ...
    *Exp_MISO_NLoS_FF_ap; 

TheoChannel_FF_ap(user_index,:)=length(UETX_dis)*Exp_MISO_FF_ap; 
end


end