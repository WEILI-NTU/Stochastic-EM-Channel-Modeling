%==========================================================================
% ==========Stochastic Green's function based channel modeling=============
%==========================================================================
%  Author: Wei Li;
%  Date:   15/05/ 2024,  NTU;
%  Version: V1.0 

%  Note: This function computes the capacity varies with T-R distance

% Please cite the paper properly if you use any part of this code


%%========================================================
clc
clear   
%%========================================================  
dis_rt_Range=[2:2:50]; % distance between the transmitter and receiver
flagLoS=1; % LoS components exist for 1
flagNLoS= 1; % NLOS components exist for 1  
%% Simulation setting
MC_max=5000; % Monto-Carlo
SNR_range=30;
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
tra.num_x=6; % number of antennas in x-direction
tra.num_y=6; % number of antennas in y-direction 
% receive antenna array setting
rec.spac=0.4*par.wav;% spacing between adjacent antennas in array
rec.num_x=4; % number of antennas in x-direction
rec.num_y=4; % number of antennas in y-direction
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
% RayleighDis=0.62*sqrt(tra.len.h^3/par.wav)/par.wav;


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

if exist(['MatrixZ_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav*10),'_',num2str(par.freq/1e9),'.txt'])==0
    error('Please run main_CST_CouplingMatrixGen.m to obtain the impedance matrix');
else
    MatrixZ=readmatrix(['MatrixZ_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav*10),'_',num2str(par.freq/1e9)]);
end

if exist(['MatrixS_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav*10),'_',num2str(par.freq/1e9),'.txt'])==0
    error('Please run main_CST_CouplingMatrixGen.m to obtain the S-parameters');
else
    MatrixS=readmatrix(['MatrixS_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav*10),'_',num2str(par.freq/1e9)]);
end

if exist(['MatrixTotalEffic_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav*10),'_',num2str(par.freq/1e9),'.txt'])==0
    error('Please run main_CST_CouplingMatrixGen.m to obtain the total efficiencies');
else
    MatrixTotalEffic=readmatrix(['MatrixTotalEffic_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav*10),'_',num2str(par.freq/1e9)]);
end

MatrixZLinear=10.^(MatrixZ./10); % dB to linear
MatrixSLinear=10.^(MatrixS./10); % dB to linear
 %% Generate mutual coupling matrix
ZS=diag(diag(MatrixZLinear));
CoupleMatrix=MatrixZLinear*inv(MatrixZLinear+ZS);


for iter_DisRT=1:length(dis_rt_Range)
    dis_rt=dis_rt_Range(iter_DisRT);  
    users.start.z=dis_rt*par.wav;


%% Distance between each pair of transmit-receiver antennas
[par,temp_UETX_dis]=DistanceCompute(par,tra,rec,users); % stores in par.UETX


% Clarke's spatial autocorrelation function  
ClarkeCor= besselj(0,2*pi*(par.UETX.distance(1,:)-par.UETX.distance(1,1) )/par.wav); 
ClarkeCor=sqrt(ClarkeCor);

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
ExpLosOper=cellfun(@(z)  norm(z,'fro').^2, CoGreenOper,'UniformOutput',false);
ExpNLosOper=cellfun(@(x,y)  (norm(termQV*x,'fro').^2+ termQV^2*norm(y,'fro')),IncoGreenExpOper,IncoGreenVarOper,'UniformOutput',false);
RatioLosNlos=cellfun(@(x,y)   (x/y),ExpLosOper,ExpNLosOper,'UniformOutput',false);

 
for snrVal=1:length(varNoise)
    tempIsoXX=0;
    tempIsoDP=0;
    tempIsoTP=0;
    tempCouXX=0;
    tempCouDP=0;
    tempCouTP=0;
    tempClarke=0;
    tempRayleigh=0;
    for MC=1:MC_max
        NLOSOper= cellfun(@(x,y)  termConstant*termQV*x+ (termConstant*termQV)*sqrt(y).*randn(size(y,1),size(y,2)),IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false);
        IsoChOper_temp= cellfun(@(x,y)  termConstant*flagLoS*x+ flagNLoS*y,CoGreenOper,NLOSOper ,'UniformOutput',false); 
        % IsoChOper_temp= cellfun(@(z) normalize(z,"norm"), IsoChOper_temp,'UniformOutput',false); 
        
        % Isolated channel capacity
        UnCoupleChOper= cellfun(@(z)z,IsoChOper_temp,'UniformOutput',false);
        UnCoupCapa=cellfun(@(x) sum(log2(1+(eig(...
            x*x'/varNoise(snrVal) ) ) ) ),  UnCoupleChOper, 'UniformOutput',false);
        tempIsoXX=tempIsoXX+UnCoupCapa{1};
        tempIsoDP=tempIsoDP+UnCoupCapa{1}+UnCoupCapa{2};
        tempIsoTP=tempIsoTP+UnCoupCapa{1}+UnCoupCapa{2}+UnCoupCapa{3}; 

        % Coupled channel (coupling matrix at transmitter, resulting in reduced energy)
        CoupleChOper= cellfun(@(z)z*CoupleMatrix*diag(MatrixTotalEffic),IsoChOper_temp,'UniformOutput',false);
        CoupCapa=cellfun(@(x) sum(log2(1+(eig(...
            x*x'/varNoise(snrVal) ) ) ) ),  CoupleChOper, 'UniformOutput',false);
        tempCouXX=tempCouXX+CoupCapa{1};
        tempCouDP=tempCouDP+CoupCapa{1}+CoupCapa{2};
        tempCouTP=tempCouTP+CoupCapa{1}+CoupCapa{2}+CoupCapa{3}; 
        
        % Rayleigh channel 
        h_rayleigh = sqrt(1/2) *(randn(rec.totalNum,tra.totalNum) + 1i*randn(rec.totalNum,tra.totalNum)); 
        h_rayleigh=h_rayleigh./(par.UETX.distance.^1);
        ClarkeChannel=h_rayleigh*ClarkeCor';
        % h_rayleigh=normalize(h_rayleigh,"norm");
        % ClarkeChannel=normalize(ClarkeChannel,"norm"); 
        tempClarke=tempClarke+sum(log2(1+(eig( ClarkeChannel*ClarkeChannel'/varNoise(snrVal) ) ) ) );
        tempRayleigh=tempRayleigh+sum(log2(1+(eig( h_rayleigh*h_rayleigh'/varNoise(snrVal) ) ) ) );  
        
        
    end
    Capa_XX_iso(iter_DisRT,snrVal)=tempIsoXX/MC_max;
    Capa_DP_iso(iter_DisRT,snrVal)=tempIsoDP/MC_max;
    Capa_TP_iso(iter_DisRT,snrVal)=tempIsoTP/MC_max; 

    Capa_XX_Cou(iter_DisRT,snrVal)=tempCouXX/MC_max;
    Capa_DP_Cou(iter_DisRT,snrVal)=tempCouDP/MC_max;
    Capa_TP_Cou(iter_DisRT,snrVal)=tempCouTP/MC_max; 

    Cap_Clarke(iter_DisRT,snrVal)=tempClarke/MC_max; 
    Cap_Rayleigh(iter_DisRT,snrVal)=tempRayleigh/MC_max;

end
 
end

figure
hold on
plot(dis_rt_Range,Capa_XX_iso(:,snrVal),'-o','color',[0.8500 0.3250 0.0980],...
    'linewidth',2,'MarkerSize',4)
plot(dis_rt_Range,Capa_XX_Cou(:,snrVal),':o','color',[0.8500 0.3250 0.0980],...
    'linewidth',2,'MarkerSize',4)
plot(dis_rt_Range,Capa_DP_iso(:,snrVal),'*-','color',[0.93,0.69,0.13],...
    'linewidth',2,'MarkerSize',4)
plot(dis_rt_Range,Capa_DP_Cou(:,snrVal),':*','color',[0.93,0.69,0.13],...
    'linewidth',2,'MarkerSize',4)
plot(dis_rt_Range,Capa_TP_iso(:,snrVal),'-d','Color',[0.13,0.62,0.77]...
    ,'linewidth',2,'MarkerSize',4) 
plot(dis_rt_Range,Capa_TP_Cou(:,snrVal),':d','Color',[0.13,0.62,0.77]...
    ,'linewidth',2,'MarkerSize',4)
plot(dis_rt_Range,Cap_Rayleigh(:,snrVal),'-^','Color',[0.42,0.73,0.36],'linewidth',2,'MarkerSize',3)
 
xlabel('distance ($\lambda$) ','Interpreter','latex');
ylabel('Capacity (bits/s/Hz)','Interpreter','latex');
legend('XX (w/o MC)','XX (with MC)','DP (w/o MC)','DP (with MC)','TP (w/o MC)',...
    'TP (with MC)','i.i.d. Rayleigh fading'...  
    ,  'Interpreter','latex');
grid on 
title([num2str([flagLoS flagNLoS]),' Capacity for ',num2str(tra.totalNum),'$\times$', num2str(rec.totalNum),' array'],'Interpreter','latex')
% xlim([dis_rt_Range(1) dis_rt_Range(end)]) 
% xticks(dis_rt_Range)

 
fprintf('end')



