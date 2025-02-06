%==========================================================================
% ==========Stochastic Green's function based channel modeling=============
%==========================================================================
%  Author: Wei Li;
%  Date:   15/05/ 2024,  NTU;
%  Version: V1.0 

%  Note: This function computes the eigenvalues of stochastic
%  Green's function based channel for differnt K-factor 

% Please cite the paper properly if you use any part of this code
%%========================================================
clc
clear   
%%======================================================== 
dis_rt_Range=[2:5:100]; % distance between the transmitter and receiver
flagLOS=1; % LoS components exist for 1
flagNLOS=0; % NLOS components exist for 1
KFactor=100; % Rician chanel 
%% Simulation setting
MC_max=1000; % Monto-Carlo
% SNR_range=-10:10;
%%  Set the parameters
par.Freq=[4,6]*1e9; % frequency range(GHz)
par.freq=(par.Freq(2)-par.Freq(1))/2+par.Freq(1); % operating frequency
par.wav=3e8/(par.freq);  % wavelength  3/(10*par.freq)?

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


% RayleighDis=2*(tra.aperture+rec.aperture)^2/par.wav; 
% RayleighDis=min(pi*tra.len.h*tra.len.v/(par.wav^2),pi*rec.len.h*rec.len.v/(par.wav^2))


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

lambda_XX_sum=zeros(length(dis_rt_Range),1);
lambda_YY_sum=zeros(length(dis_rt_Range),1);
lambda_ZZ_sum=zeros(length(dis_rt_Range),1);
lambda_XY_sum=zeros(length(dis_rt_Range),1);
lambda_XZ_sum=zeros(length(dis_rt_Range),1);
lambda_YZ_sum=zeros(length(dis_rt_Range),1);
lambda_Rician_sum=zeros(length(dis_rt_Range),1);
for iter_DisRT=1:length(dis_rt_Range)
    dis_rt=dis_rt_Range(iter_DisRT);  
    users.start.z=dis_rt*par.wav;

%% Distance between each pair of transmit-receiver antennas
[par,temp_UETX_dis]=DistanceCompute(par,tra,rec,users); % stores in par.UETX

%% Isolated channle modeling
%% Cavity spectrum: \lambda_i
% alpha=par.wavenumber^3*par.cavity/(2*pi^2*par.Qfactor);
% [cavity_ki,CavityEigTerm]=GOE_Eigen(num_M,alpha,par);

% Dyadic Green's function and stochastic Green's function
 % Incoherent part 
% [IncoGreenExp,IncoGreenVar]=IncoGreenFuncIndi(par,tra,rec,cavity_ki);
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
 
    tempIsoXX=0;
    tempIsoDP=0;
    tempIsoTP=0;
    for MC=1:MC_max
        NLOSOper= cellfun(@(x,y)  termConstant*x+ (termConstant)*sqrt(y).*randn(size(y,1),size(y,2)),IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false);
        IsoChOper= cellfun(@(x,y,fl,nl)  termConstant*flagLOS*x+ flagNLOS*y,CoGreenOper,NLOSOper ,'UniformOutput',false);
        % normalization
        % IsoChOper= cellfun(@(z) normalize(z,"norm"), IsoChOper,'UniformOutput',false);
        SumIsoChanOper= cellfun(@(z) z*z', IsoChOper,'UniformOutput',false); 
        % [~,lambda_XX_temp,~]=svd(SumIsoChanOper{1});
        % [~,lambda_DP_temp,~]=svd(SumIsoChanOper{1}+SumIsoChanOper{2});
        % [~,lambda_TP_temp,~]=svd(SumIsoChanOper{1}+SumIsoChanOper{2}+SumIsoChanOper{3});
        lambda_XX=sort(eig(SumIsoChanOper{1}),'descend');
        lambda_YY=sort(eig(SumIsoChanOper{2}),'descend');
        lambda_ZZ=sort(eig(SumIsoChanOper{3}),'descend');
        lambda_XY=sort(eig(SumIsoChanOper{4}),'descend');
        lambda_XZ=sort(eig(SumIsoChanOper{5}),'descend');
        lambda_YZ=sort(eig(SumIsoChanOper{6}),'descend');

        lambda_XX_sum(iter_DisRT,1)=lambda_XX_sum(iter_DisRT,1)+sum(lambda_XX);
        lambda_YY_sum(iter_DisRT,1)=lambda_YY_sum(iter_DisRT,1)+sum(lambda_YY);
        lambda_ZZ_sum(iter_DisRT,1)=lambda_ZZ_sum(iter_DisRT,1)+sum(lambda_ZZ);
        lambda_XY_sum(iter_DisRT,1)=lambda_XY_sum(iter_DisRT,1)+sum(lambda_XY);
        lambda_XZ_sum(iter_DisRT,1)=lambda_XZ_sum(iter_DisRT,1)+sum(lambda_XZ);
        lambda_YZ_sum(iter_DisRT,1)=lambda_YZ_sum(iter_DisRT,1)+sum(lambda_YZ);

        %% Rician channel
        LoSChannel= cellfun(@(x)  termConstant*x,CoGreenOper,'UniformOutput',false);
        % LoSChannel= cellfun(@(z) normalize(z,"norm"), LoSChannel,'UniformOutput',false);
        SumLoSChanOper= cellfun(@(z) z*z', LoSChannel,'UniformOutput',false); 
        
        Corician=real((1-1/(KFactor+1))*SumLoSChanOper{1}+ ...
             1/(KFactor+1)*exp(1j*2*pi*dis_rt_Range(iter_DisRT) ));
        lambda_RicianTemp=sort(eig(Corician),'descend');
        lambda_Rician_sum(iter_DisRT,1)=lambda_Rician_sum(iter_DisRT,1)+sum(lambda_RicianTemp);
    end 
end
lambda_XX_MC= (lambda_XX_sum)/MC_max;
lambda_YY_MC= (lambda_YY_sum)/MC_max;
lambda_ZZ_MC= (lambda_ZZ_sum)/MC_max;
lambda_XY_MC= (lambda_XY_sum)/MC_max;
lambda_XZ_MC= (lambda_XZ_sum)/MC_max;
lambda_YZ_MC= (lambda_YZ_sum)/MC_max;
lambda_Rician_MC=lambda_Rician_sum/MC_max;

figure 
plot(dis_rt_Range,10*log10(lambda_XX_MC),'-*','Color',[0.86,0.27,0.27],'linewidth',2)
hold on
plot(dis_rt_Range,10*log10(lambda_YY_MC),'--o','Color',[0.35,0.35,0.87],'linewidth',2)
plot(dis_rt_Range,10*log10(lambda_ZZ_MC),'o-','Color',[0.92,0.42,0.92],'linewidth',2)
plot(dis_rt_Range,10*log10(lambda_XY_MC),'-^','color',[0.27,0.73,0.87],'linewidth',2)
plot(dis_rt_Range,10*log10(lambda_XZ_MC),'-v','color',[0.87,0.66,0.16],'linewidth',2)
plot(dis_rt_Range,10*log10(lambda_YZ_MC),'-d','Color',[0.78,0.58,0.87],'linewidth',2)
% plot(dis_rt_Range,10*log10(lambda_Rician_MC),'k--','linewidth',2)
xlabel('distance ($\lambda$)','Interpreter','latex');
ylabel('Eigenvalue (dB)','Interpreter','latex');
legend({'XX','YY','ZZ','XY','XZ','YZ' } ,'Interpreter','latex','Location','NorthEast');
grid on 
title([ num2str([flagLOS flagNLOS]),', $N_s$=',num2str(tra.totalNum),', $N_r$=',num2str(rec.totalNum), ...
    ', $\Delta_s$=',num2str(tra.spac/par.wav),'$\lambda$', ...
    ', $\Delta_r$=',num2str(rec.spac/par.wav),'$\lambda$'],'Interpreter','latex')
 
 

% xline(RayleighDis/par.wav,'--',{'Rayleigh','distance'},...
%     'LabelOrientation','horizontal','LabelHorizontalAlignment','right',...
%     'LabelVerticalAlignment','middle',...
%     'linewidth',2,'Interpreter','latex','HandleVisibility','off')
% 
% 
 



fprintf('end')



