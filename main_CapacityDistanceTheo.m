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
dis_rt_Range=[2:10:300]; % distance between the transmitter and receiver
flagLoS=1; % LoS components exist for 1
flagNLoS=1; % NLOS components exist for 1   
%% Simulation setting
MC_max=5000; % Monto-Carlo
SNR_range= -10;
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
tra.spac=0.4*par.wav;% spacing between adjacent antennas in array 
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
 

rec.len.h=rec.num_x*par.rad*2+(rec.num_x-1)*rec.spac;
rec.len.v=rec.num_y*par.len+(rec.num_y-1)*rec.spac;
rec.aperture=sqrt(rec.len.h^2+rec.len.v^2);
rec.totalNum=rec.num_x*rec.num_y;
rec.indi.x=par.rad*2+rec.spac;
rec.indi.y=par.len+rec.spac; 


% RayleighDis=2*(tra.aperture)^2/par.wav;  
 
 RayleighDis=0.62*sqrt(tra.len.h^3/par.wav)/par.wav;
 

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

  
for snrVal=1:length(varNoise)
    tempIsoXX=0;
    tempIsoDP=0;
    tempIsoTP=0;
    tempCouXX=0;
    tempCouDP=0;
    tempCouTP=0; 
    for MC=1:MC_max
        NLOSOper= cellfun(@(x,y)  termConstant*termQV*x+ (termConstant*termQV)*sqrt(y).*randn(size(y,1),size(y,2)),IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false);
        IsoChOper_temp= cellfun(@(x,y)  termConstant*flagLoS*x+ flagNLoS*y,CoGreenOper,NLOSOper ,'UniformOutput',false); 
     
        % Isolated channel capacity
        UnCoupleChOper= cellfun(@(z)z,IsoChOper_temp,'UniformOutput',false);
        UnCoupCapa=cellfun(@(x) sum(log2(1+(eig(...
            x*x'/varNoise(snrVal) ) ) ) ),  UnCoupleChOper, 'UniformOutput',false);
        tempIsoXX=tempIsoXX+UnCoupCapa{1};
        tempIsoDP=tempIsoDP+UnCoupCapa{1}+UnCoupCapa{2};
        tempIsoTP=tempIsoTP+UnCoupCapa{1}+UnCoupCapa{2}+UnCoupCapa{3}; 

 
 
    end
    Capa_XX_iso(iter_DisRT,snrVal)=tempIsoXX/MC_max;
    Capa_DP_iso(iter_DisRT,snrVal)=tempIsoDP/MC_max;
    Capa_TP_iso(iter_DisRT,snrVal)=tempIsoTP/MC_max; 
 
end
 
[TheoChannel_Ac,TheoChannel_FF,TheoChannel_FF_ap]=ChApproximate(par,temp_UETX_dis,termConstant,flagLoS,flagNLoS);
Capa_Theo_Ac(iter_DisRT,snrVal)=sum(log2(1+((sum(TheoChannel_Ac)/varNoise(snrVal) ) ) ) ); 
Capa_Theo_FF(iter_DisRT,snrVal)=sum(log2(1+((sum(TheoChannel_FF)/varNoise(snrVal) ) ) ) ); 
Capa_Theo_FF_ap(iter_DisRT,snrVal)=sum(log2(1+((sum(TheoChannel_FF_ap)/varNoise(snrVal) ) ) ) ); 

end

figure
hold on
plot(dis_rt_Range,Capa_XX_iso(:,snrVal),'-*','color',[0.91,0.47,0.47],...
    'linewidth',2,'MarkerSize',4) 
% plot(dis_rt_Range,Capa_Theo_Ac(:,snrVal),'-x','Color',[0.47,0.67,0.19]...
%     ,'linewidth',2,'MarkerSize',1)
plot(dis_rt_Range,Capa_Theo_FF(:,snrVal),'-o','Color',[0.21,0.69,0.62]...
    ,'linewidth',2,'MarkerSize',5)
plot(dis_rt_Range,Capa_Theo_FF_ap(:,snrVal),'-d','Color',[0.45,0.43,0.86]...
    ,'linewidth',2,'MarkerSize',5)
 
xlabel('distance ($\lambda$) ','Interpreter','latex');
ylabel('Capacity (bits/s/Hz)','Interpreter','latex');
% legend('XX (Ideal)','MISO (Ideal, Accurate)','MISO (Ideal, Far field)',...  
%     'MISO (Ideal, Far field, small aperture)' ,'Interpreter','latex');
legend('Practical', 'Theo. (lower-bound)',...  
    'Theo. (upper-bound)' ,'Interpreter','latex');
grid on 
title([num2str([flagLoS flagNLoS]),' Capacity for ',num2str(tra.totalNum),'$\times$', num2str(rec.totalNum), ...
    ' array with spacing ',num2str(tra.spac/par.wav),'$\lambda$'],'Interpreter','latex')
% xlim([dis_rt_Range(1) dis_rt_Range(end)]) 
% xticks(dis_rt_Range)

xline(RayleighDis,'--', ...
    'LabelOrientation','horizontal','LabelHorizontalAlignment','left',...
    'LabelVerticalAlignment','middle',...
    'linewidth',2,'Interpreter','latex','HandleVisibility','off')
 

fprintf('end')



