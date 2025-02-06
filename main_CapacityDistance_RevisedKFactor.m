%==========================================================================
% ==========Stochastic Green's function based channel modeling=============
%==========================================================================
%  Author: Wei Li;
%  Date:   03/07/ 2024,  NTU;
%  Version: V1.0 

%  Note: This function computes the capacity varies with T-R distance

% Please cite the paper properly if you use any part of this code


%%========================================================
clc
clear   
%%========================================================  
dis_rt_Range=[2:5:100]; % distance between the transmitter and receiver
 
KFactor=10; % K factor for Rician channel 
%% Simulation setting
MC_max=5000; % Monto-Carlo
SNR_range=-10;
%%  Set the parameters
par.Freq=[4,6]*1e9; % frequency range(GHz)
par.freq=(par.Freq(2)-par.Freq(1))/2+par.Freq(1); % operating frequency
par.wav=3e8/(par.freq); % wavelength  3/(10*par.freq)?
par.wavenumber=2*pi/par.wav;
 

alpha= 1e3;
% par.cavity=2*pi^2*par.Qfactor*alpha/(par.wavenumber^3);
par.cavity=(400*par.wav)^3; 
par.Qfactor= par.wavenumber^3*par.cavity/(2*pi^2*alpha); 
num_M=2000; % number of cavity eigenvalues
par.mu=4*pi*1e-7;   
par.omega=2*pi*par.freq;
par.curr=1; 
% dipole setting
par.rad=0.01*par.wav; % radius of dipole
par.len=0.43*par.wav; % length of dipole % dipole size:radius*len
par.gap=0.01*par.wav; % The gap inside the dipole
% transmit antenna array setting
tra.spac=0.1*par.wav;% spacing between adjacent antennas in array 
tra.num_x=100; % number of antennas in x-direction
tra.num_y=1; % number of antennas in y-direction 
% receive antenna array setting
rec.spac=0.4*par.wav;% spacing between adjacent antennas in array
rec.num_x=1; % number of antennas in x-direction
rec.num_y=1; % number of antennas in y-direction
rec.total=rec.num_x*rec.num_y; % total number of antennas in HMIMO
% transmission setting
% HMIMOS setting   
tra.len.h=tra.num_x*par.rad*2+(tra.num_x-1)*tra.spac;
tra.len.v=tra.num_y*par.len+(tra.num_y-1)*tra.spac;
tra.aperture=sqrt(tra.len.h^2+tra.len.v^2);
tra.totalNum=tra.num_x*tra.num_y;
tra.indi.x=par.rad*2+tra.spac;
tra.indi.y=par.len+tra.spac;
 

rec.len.h=rec.num_x*par.rad*2+(rec.num_x-1)*tra.spac;
rec.len.v=rec.num_y*par.len+(rec.num_y-1)*tra.spac;
rec.aperture=sqrt(rec.len.h^2+rec.len.v^2);
rec.totalNum=rec.num_x*rec.num_y;
rec.indi.x=par.rad*2+rec.spac;
rec.indi.y=par.len+rec.spac; 


RayleighDis=2*(tra.aperture+rec.aperture)^2/par.wav;  

 

tra.start.h=0;
tra.start.v=0;
tra.coordinateX=tra.start.h:tra.indi.x...
        :tra.start.h+(tra.num_x-1)*tra.indi.x;
tra.coordinateX=kron(ones(1,tra.num_y),tra.coordinateX);
tra.coordinateY=tra.start.v:tra.indi.y...
        :tra.start.v+(tra.num_y-1)*tra.indi.y;
tra.coordinateY=kron(tra.coordinateY,ones(1,tra.num_x));

users.num=1; 
users.start.h=0*par.wav:rec.aperture:(users.num-1)*rec.aperture; 
users.start.v=0:rec.aperture:(users.num-1)*rec.aperture;  

varNoise=10.^(-0.1.*SNR_range);
PowerCurrent=eye(tra.totalNum);

 
for iter_DisRT=1:length(dis_rt_Range)
    dis_rt=dis_rt_Range(iter_DisRT);  
    users.start.z=dis_rt*par.wav;


%% Distance between each pair of transmit-receiver antennas
[par,temp_UETX_dis]=DistanceCompute(par,tra,rec,users); % stores in par.UETX

%% Isolated channle modeling

% Dyadic Green's function and stochastic Green's function
 % Incoherent part  
[IncoGreenExp,IncoGreenVar]=IncoGreenFunc(par,tra,rec,users);
IncoGreenExpOper=cellfun(@(z) z,struct2cell(IncoGreenExp) ,'UniformOutput',false);
IncoGreenVarOper=cellfun(@(z) z,struct2cell(IncoGreenVar) ,'UniformOutput',false);
% Coherent channel 
[CoGreen]=CoGreenFunc(par,temp_UETX_dis);
CoGreenOper=cellfun(@(z) real(z),struct2cell(CoGreen) ,'UniformOutput',false); 
 
% termQV=CavityEigTerm; % Finite number of eigenvalues
termQV=pi*par.Qfactor/(par.wavenumber^2); % Infinite waves
% termConstant= 1i*2*pi*par.freq*4*pi*10e-7;
termConstant=1i*par.wavenumber*120*pi;

% Compute the ratio between the LOS and NLOS 
ExpLosOper=cellfun(@(z)  norm(termConstant*z,'fro').^2, CoGreenOper,'UniformOutput',false);
ExpNLosOper=cellfun(@(x,y)  (norm(termConstant*x,'fro').^2+norm((termConstant)^2*y,'fro')),IncoGreenExpOper,IncoGreenVarOper,'UniformOutput',false);
RatioLosNlos=cellfun(@(x,y)   (x/y),ExpLosOper,ExpNLosOper,'UniformOutput',false);
K_constant=cellfun(@(z)  KFactor/z, RatioLosNlos,'UniformOutput',false); 
flagNLoS=cellfun(@(z) sqrt(1/(1+z)), K_constant,'UniformOutput',false);  
flagLoS=cellfun(@(z) sqrt(1-1/(1+z)), K_constant,'UniformOutput',false); 


for snrVal=1:length(varNoise)
    tempIsoXX=0;
    tempIsoDP=0;
    tempIsoTP=0;
    tempCouXX=0;
    tempCouDP=0;
    tempCouTP=0;
    Cap_Rician_sum=0;
    for MC=1:MC_max
        NLOSOper= cellfun(@(x,y)  termConstant*termQV*x+ (termConstant*termQV)*sqrt(y).*randn(size(y,1),size(y,2)),IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false);
        % K-factor
        NLOSOper= cellfun(@(x,y)  x*y,NLOSOper,flagNLoS ,'UniformOutput',false);
        LOSOper= cellfun(@(x,y)  termConstant*x*y,CoGreenOper,flagLoS ,'UniformOutput',false);
        IsoChOper_temp= cellfun(@(x,y)  x+ y,LOSOper,NLOSOper, 'UniformOutput',false); 
        IsoChOper= cellfun(@(z) normalize(z,"norm"), IsoChOper_temp,'UniformOutput',false); 
        
        
        % Isolated channel capacity
        UnCoupleChOper= cellfun(@(z)z,IsoChOper_temp,'UniformOutput',false);
        UnCoupCapa=cellfun(@(x) sum(log2(1+(eig(...
            x*x'/varNoise(snrVal) ) ) ) ),  UnCoupleChOper, 'UniformOutput',false);
        tempIsoXX=tempIsoXX+UnCoupCapa{1};
        tempIsoDP=tempIsoDP+UnCoupCapa{1}+UnCoupCapa{2};
        tempIsoTP=tempIsoTP+UnCoupCapa{1}+UnCoupCapa{2}+UnCoupCapa{3}; 


         %% Rician channel
        LoSChannel= cellfun(@(x)  termConstant*x,CoGreenOper,'UniformOutput',false);
        % Compute the ratio between the LOS and NLOS  
        flagNLoS_Ric=  sqrt(1/(1+KFactor)) ; 
        flagLoS_Ric= sqrt(1-1/(1+KFactor)) ; 
        % To generate a Rayleigh fading channel gain
        h_rayleigh = sqrt(1/2) *(randn(rec.num_x*rec.num_y, ...
            tra.num_x*tra.num_y) + 1i*randn(rec.num_x*rec.num_y,tra.num_x*tra.num_y)); 
        % To generate a Rician fading channel gain
         h_rician = flagLoS_Ric*flagLoS{2} *  (LoSChannel{2}) +flagNLoS_Ric*  h_rayleigh; 
        Cap_Rician_sum=Cap_Rician_sum+sum(log2(1+(eig(h_rician*h_rician'/varNoise(snrVal) ) ) ) );
 
    end
    Capa_XX_iso(iter_DisRT,snrVal)=tempIsoXX/MC_max;
    Capa_DP_iso(iter_DisRT,snrVal)=tempIsoDP/MC_max;
    Capa_TP_iso(iter_DisRT,snrVal)=tempIsoTP/MC_max; 

    Capa_Rician(iter_DisRT,snrVal)=Cap_Rician_sum/MC_max;
end

% Theoretical channel capacity of MISO systems 
% Theo_CH_all= cellfun(@(x,y,z)  termConstant*flagLoS*x+ termConstant*termQV*y+(termConstant*termQV)*sqrt(z),CoGreenOper,IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false); 
% Theo_Ch_cap=cellfun(@(x) sum(log2(1+(eig(...
%         x*x'/varNoise(snrVal) ) ) ) ),  Theo_CH_all, 'UniformOutput',false);

% Theo_CH_all= cellfun(@(x,y,z)  flagLoS*norm(termConstant)^2*x.*conj(x) ...
%     + flagNLoS*norm(termConstant*termQV)^2*y.*conj(y) ... 
%     +flagNLoS*norm(termConstant*termQV)*z,CoGreenOper,IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false); 
% Theo_Ch_cap= cellfun(@(z)sum(z),Theo_CH_all,'UniformOutput',false);
% Theo_Ch_cap=cellfun(@(x) sum(log2(1+(eig(...
%         x/varNoise(snrVal) ) ) ) ),  Theo_Ch_cap, 'UniformOutput',false);
% 
% Capa_Theo(iter_DisRT,snrVal)=Theo_Ch_cap{1};

[TheoChannel_Ac,TheoChannel_FF,TheoChannel_FF_ap]=ChApproximate(par,temp_UETX_dis,termConstant,flagLoS{1},flagNLoS{1});
Capa_Theo_Ac(iter_DisRT,snrVal)=sum(log2(1+(eig(sum(TheoChannel_Ac)/varNoise(snrVal) ) ) ) ); 
Capa_Theo_FF(iter_DisRT,snrVal)=sum(log2(1+(eig(sum(TheoChannel_FF)/varNoise(snrVal) ) ) ) ); 
Capa_Theo_FF_ap(iter_DisRT,snrVal)=sum(log2(1+(eig(sum(TheoChannel_FF_ap)/varNoise(snrVal) ) ) ) ); 

end

figure
hold on
plot(dis_rt_Range,Capa_XX_iso(:,snrVal),'-','color',[0.13,0.62,0.77],...
    'linewidth',2,'MarkerSize',1) 
plot(dis_rt_Range,Capa_Theo_Ac(:,snrVal),'-x','Color',[0.47,0.67,0.19]...
    ,'linewidth',2,'MarkerSize',1)
plot(dis_rt_Range,Capa_Theo_FF(:,snrVal),'-x','Color',[0.90,0.61,0.45]...
    ,'linewidth',2,'MarkerSize',1)
plot(dis_rt_Range,Capa_Theo_FF_ap(:,snrVal),'-x','Color',[0.56,0.44,0.93]...
    ,'linewidth',2,'MarkerSize',1)
plot(dis_rt_Range,Capa_Rician(:,snrVal),'-kx'...
    ,'linewidth',2,'MarkerSize',1)
xlabel('distance ($\lambda$) ','Interpreter','latex');
ylabel('Capacity (bits/s/Hz)','Interpreter','latex');
legend('XX (Ideal)','MISO (Ideal, Accurate)','MISO (Ideal, Far field)',...  
    'MISO (Ideal, Far field, small aperture)',['Rician Channel (K=',num2str(KFactor),')'],'Interpreter','latex');
grid on 
title(['$K=$',num2str(KFactor),' Capacity for ',num2str(tra.totalNum),'$\times$', num2str(rec.totalNum),' array'],'Interpreter','latex')
% xlim([dis_rt_Range(1) dis_rt_Range(end)]) 
% xticks(dis_rt_Range)

% xline(RayleighDis/par.wav,'--',{'Rayleigh','distance'},...
%     'LabelOrientation','horizontal','LabelHorizontalAlignment','right',...
%     'LabelVerticalAlignment','middle',...
%     'linewidth',2,'Interpreter','latex','HandleVisibility','off')
 

fprintf('end')



