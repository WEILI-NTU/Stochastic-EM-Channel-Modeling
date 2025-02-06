%==========================================================================
% ==========Stochastic Green's function based channel modeling=============
%==========================================================================
%  Author: Wei Li;
%  Date:   15/05/ 2024,  NTU;
%  Version: V1.0 

%  Note: This function computes the autocorrelation function at 
%        the transmitter for different K values and Rician channel
% The varying distance between T-R results in varying correlation with
% spacing

% Please cite the paper properly if you use any part of this code


%%========================================================
clc
clear   
%%======================================================== 
dis_rt_Range=[3]; % distance between the transmitter and receiver
flagNLoS=0*ones(6,1);
flagLoS=1*ones(6,1);
flagNLoS=num2cell(flagNLoS);
flagLoS=num2cell(flagLoS);
KFactor=0;
 MC_max=1; % Monto-Carlo
%%  Set the parameters
par.Freq=[4,6]*1e9; % frequency range(GHz)
par.freq=(par.Freq(2)-par.Freq(1))/2+par.Freq(1); % operating frequency
par.wav=3e8/(par.freq); % wavelength  3/(10*par.freq)?
par.wavenumber=2*pi/par.wav;
alpha= 1e3;
% par.cavity=2*pi^2*par.Qfactor*alpha/(par.wavenumber^3);
par.cavity=(400*par.wav)^3; 
par.Qfactor= par.wavenumber^3*par.cavity/(2*pi^2*alpha); 
num_M=1000; % number of cavity eigenvalues
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
rec.num_x=60; % number of antennas in x-direction
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


lambda_XX_sum=zeros(( tra.totalNum),1);
lambda_YY_sum=zeros(( tra.totalNum),1);
lambda_ZZ_sum=zeros(( tra.totalNum),1);
lambda_XY_sum=zeros(( tra.totalNum),1);
lambda_XZ_sum=zeros(( tra.totalNum),1);
lambda_YZ_sum=zeros(( tra.totalNum),1);
lambda_Rician_sum=zeros(( tra.totalNum),1);
for iter_DisRT=1:length(dis_rt_Range)
    dis_rt=dis_rt_Range(iter_DisRT);  
    users.start.z=dis_rt*par.wav;


%% Distance between each pair of transmit-receiver antennas
[par,temp_UETX_dis]=DistanceCompute(par,tra,rec,users); % stores in par.UETX


% Clarke's spatial autocorrelation function  
ClarkeCor= besselj(0,2*pi*(par.UETX.distance(1,:)-par.UETX.distance(1,1) )/par.wav); 


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
termConstant= 1i*2*pi*par.freq*4*pi*10e-7;

%  K_constant= cellfun(@(x,y,z)   ((norm(termConstant*x*x'))) ...
%     /((norm(termConstant*termQV*y*y')) ... 
%     +(norm(termConstant*termQV*z)))  ...
%     ,CoGreenOper,IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false); 
% flagNLoS= cellfun(@(x) sqrt(x/(KFactor+x)) ,K_constant,'UniformOutput',false); 
% flagLoS= cellfun(@(x) sqrt(1-x^2) ,flagNLoS,'UniformOutput',false);  
  

 

Cor_Rician_sum=0;
for MC=1:MC_max
    % K-factor
    NLOSOper= cellfun(@(x,y)  termConstant*termQV*x+   (termConstant*termQV)*sqrt(y).*randn(size(y,1),size(y,2)),IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false);
    NLOSOper= cellfun(@(x,y)  x*y,NLOSOper,flagNLoS ,'UniformOutput',false);
    LOSOper= cellfun(@(x)  termConstant*x,CoGreenOper ,'UniformOutput',false);
    LOSOper= cellfun(@(x,y)  x*y,LOSOper,flagLoS ,'UniformOutput',false);
    IsoChOper_temp= cellfun(@(x,y) x+ y,LOSOper,NLOSOper, 'UniformOutput',false); 
    % IsoChOper_temp= cellfun(@(x,y,z) sqrt(KFactor/z)*x+  y,LOSOper,NLOSOper,K_constant, 'UniformOutput',false); 

    % Isolated channel capacity
    UnCoupleChOper= cellfun(@(z)z,IsoChOper_temp,'UniformOutput',false);
    % UnCoupleChOper= cellfun(@(z) normalize(z,"norm"), IsoChOper_temp,'UniformOutput',false); 
    UnCoupCapa=cellfun(@(x)  x'*x,  UnCoupleChOper, 'UniformOutput',false); 
     
     lambda_XX=sort(abs(eig(UnCoupCapa{1})),'descend');
        lambda_YY=sort(abs(eig(UnCoupCapa{2})),'descend');
        lambda_ZZ=sort(abs(eig(UnCoupCapa{3})),'descend');
        lambda_XY=sort(abs(eig(UnCoupCapa{4})),'descend');
        lambda_XZ=sort(abs(eig(UnCoupCapa{5})),'descend');
        lambda_YZ=sort(abs(eig(UnCoupCapa{6})),'descend');

        lambda_XX_sum =lambda_XX_sum + (lambda_XX);
        lambda_YY_sum =lambda_YY_sum + (lambda_YY);
        lambda_ZZ_sum =lambda_ZZ_sum + (lambda_ZZ);
        lambda_XY_sum =lambda_XY_sum + (lambda_XY);
        lambda_XZ_sum =lambda_XZ_sum + (lambda_XZ);
        lambda_YZ_sum =lambda_YZ_sum + (lambda_YZ);

        % %% Rician channel
        % LoSChannel= cellfun(@(x)  termConstant*x,CoGreenOper,'UniformOutput',false);
        % % LoSChannel= cellfun(@(z) normalize(z,"norm"), LoSChannel,'UniformOutput',false);
        % SumLoSChanOper= cellfun(@(z) z*z', LoSChannel,'UniformOutput',false); 
        % KFactor=10;
        % Corician=real(SumLoSChanOper{1}+KFactor ...
        %     *exp(1j*2*pi*dis_rt_Range(iter_DisRT)*par.wav/20))/(KFactor+1);
        % lambda_RicianTemp=sort(eig(Corician),'descend');
        % lambda_Rician_sum(iter_DisRT,1)=lambda_Rician_sum(iter_DisRT,1)+sum(lambda_RicianTemp);
     %% Rician channel
        LoSChannel= cellfun(@(x)  termConstant*x,CoGreenOper,'UniformOutput',false);
        % LoSChannel= cellfun(@(z) normalize(z,"norm"), LoSChannel,'UniformOutput',false);
        SumLoSChanOper= cellfun(@(z) z'*z, LoSChannel,'UniformOutput',false); 
        
        Corician=real((1-1/(KFactor+1))*SumLoSChanOper{1}+ ...
             1/(KFactor+1)*exp(1j*2*pi*dis_rt_Range(iter_DisRT) )); 
    lambda_RicianTemp=sort(abs(eig(Corician)),'descend');
    lambda_Rician_sum =lambda_Rician_sum + (lambda_RicianTemp);
 
         
end
 
end     
lambda_XX_MC= (lambda_XX_sum)/MC_max;
lambda_YY_MC= (lambda_YY_sum)/MC_max;
lambda_ZZ_MC= (lambda_ZZ_sum)/MC_max;
lambda_XY_MC= (lambda_XY_sum)/MC_max;
lambda_XZ_MC= (lambda_XZ_sum)/MC_max;
lambda_YZ_MC= (lambda_YZ_sum)/MC_max;
lambda_Rician_MC=lambda_Rician_sum/MC_max;

XRange=par.UETX.distance(1,:)/par.wav-par.UETX.distance(1,1)/par.wav;
% XLimRange=[0 5];

% Plot figures
figure 
plot(1:tra.totalNum, 10*log10(lambda_XX_MC),'r-','linewidth',2)
hold on
plot(1:tra.totalNum, 10*log10(lambda_YY_MC),'b-','linewidth',2)
plot(1:tra.totalNum, 10*log10(lambda_ZZ_MC),'m-','linewidth',2)
plot(1:tra.totalNum, 10*log10(lambda_XY_MC),'color',[0.13,0.62,0.77],'linewidth',2)
plot(1:tra.totalNum, 10*log10(lambda_XZ_MC),'color',[0.93,0.69,0.13],'linewidth',2)
plot(1:tra.totalNum, 10*log10(lambda_YZ_MC),'Color',[0.64,0.08,0.18],'linewidth',2)
plot(1:tra.totalNum, 10*log10(lambda_Rician_MC),'k:','linewidth',2)
xlabel('Eigenvalue number ($\lambda$)','Interpreter','latex');
ylabel('Eigenvalue','Interpreter','latex');
legend({'XX','YY','ZZ','XY','XZ','YZ',['Rician ($K=$',num2str(KFactor),')']} ,'Interpreter','latex','Location','NorthEast');
grid on 
title(['[', num2str([flagLoS{1},flagNLoS{1}]),'], $N_s$=',num2str(tra.totalNum),', $N_r$=',num2str(rec.totalNum), ...
    ', $\Delta_s$=',num2str(tra.spac/par.wav),'$\lambda$', ...
    ', $\Delta_r$=',num2str(rec.spac/par.wav),'$\lambda$'],'Interpreter','latex')
 ylim([-100 150]);
 
pi*min(rec.len.h*rec.len.v,tra.len.h*tra.len.v)/(par.wav^2) 
fprintf('end')



