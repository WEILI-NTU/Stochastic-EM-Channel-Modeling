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
dis_rt=10; % distance between the transmitter and receiver
KFactor=2; % K factor for Rician channel 
MC_max=1000; % Monto-Carlo
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
tra.spac=0.05*par.wav;% spacing between adjacent antennas in array
num_Sampling=100;
tra.num_x=num_Sampling; % number of antennas in x-direction
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
users.start.h=-dis_rt*par.wav:rec.aperture:(users.num-1)*rec.aperture; 
users.start.v=0:rec.aperture:(users.num-1)*rec.aperture;  
users.start.z=0*par.wav;


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
termConstant= 1i*2*pi*par.freq*4*pi*10e-7;

%  K_constant= cellfun(@(x,y,z)   (sum(norm(termConstant)^2*x.*conj(x))) ...
%     /(sum(norm(termConstant*termQV)^2*y.*conj(y)) ... 
%     +sum(norm(termConstant*termQV)*z))/length(z) ...
%     ,CoGreenOper,IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false); 
% flagNLoS= cellfun(@(x) sqrt(x/(KFactor+x)) ,K_constant,'UniformOutput',false); 
% flagLoS= cellfun(@(x) sqrt(1-x^2) ,flagNLoS,'UniformOutput',false);  
 
% Compute the ratio between the LOS and NLOS 
ExpLosOper=cellfun(@(z)  norm(termConstant*z,'fro').^2, CoGreenOper,'UniformOutput',false);
ExpNLosOper=cellfun(@(x,y)  (norm(termConstant*x,'fro').^2+norm((termConstant)^2*y,'fro')),IncoGreenExpOper,IncoGreenVarOper,'UniformOutput',false);
RatioLosNlos=cellfun(@(x,y)   (x/y),ExpLosOper,ExpNLosOper,'UniformOutput',false);
K_constant=cellfun(@(z)  KFactor/z, RatioLosNlos,'UniformOutput',false); 
flagNLoS=cellfun(@(z) sqrt(1/(1+z)), K_constant,'UniformOutput',false);  
flagLoS=cellfun(@(z) sqrt(1-1/(1+z)), K_constant,'UniformOutput',false); 


% Clarke's spatial autocorrelation function  
ClarkeCor= besselj(0,2*pi*(par.UETX.distance(1,:)-par.UETX.distance(1,1) )/par.wav); 
% ClarkeCor=  sinc((par.UETX.distance(1,:)-par.UETX.distance(1,1) )/par.wav); 
 

 
% Generate isolated channel matrix 
tempIsoXX=0;
tempIsoDP=0;

Cor_Rician_sum=0;
for MC=1:MC_max
    % K-factor
    NLOSOper= cellfun(@(x,y)  termConstant*x+   (termConstant)*sqrt(y).*randn(size(y,1),size(y,2)),IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false);
    NLOSOper= cellfun(@(x,y)  x*y,NLOSOper,flagNLoS ,'UniformOutput',false);
    LOSOper= cellfun(@(x)  termConstant*x,CoGreenOper ,'UniformOutput',false);
    LOSOper= cellfun(@(x,y)  x*y,LOSOper,flagLoS ,'UniformOutput',false);
    IsoChOper_temp= cellfun(@(x,y) x+ y,LOSOper,NLOSOper, 'UniformOutput',false); 
    % IsoChOper_temp= cellfun(@(x,y,z) sqrt(KFactor/z)*x+  y,LOSOper,NLOSOper,K_constant, 'UniformOutput',false); 
    % IsoChOper_temp= cellfun(@(z) normalize(z,"norm"), IsoChOper_temp,'UniformOutput',false); 

    % Isolated channel capacity
    UnCoupleChOper= cellfun(@(z)z,IsoChOper_temp,'UniformOutput',false);
    UnCoupCapa=cellfun(@(x)  x'*x,  UnCoupleChOper, 'UniformOutput',false);
    tempIsoXX=tempIsoXX+UnCoupCapa{2}(1,:)/max(UnCoupCapa{2}(1,:));
    tempIsoDP=tempIsoDP+(UnCoupCapa{1}(1,:)/max(UnCoupCapa{1}(1,:)) ...
        +UnCoupCapa{2}(1,:)/max(UnCoupCapa{2}(1,:)))/2;  
     

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
    % h_rician=normalize(h_rician,"norm");
   
    Rician_temp=h_rician'*h_rician;
    Cor_Rician_sum=Cor_Rician_sum+real(Rician_temp(1,:)/max(Rician_temp(1,:)));
    

    % correlation_rician=(besselj(0,par.UETX.distance*2*pi/num_Sampling)+KFactor*cos(par.UETX.distance*2*pi/num_Sampling)+1j*KFactor*...
    %         sin(par.UETX.distance*2*pi/num_Sampling))/(1+KFactor);
    % correlation_rician=(besselj(0,par.UETX.distance*2*pi)+KFactor*cos(par.UETX.distance*2*pi)+1j*KFactor*...
    %         sin(par.UETX.distance*2*pi))/(1+KFactor);
    % Cor_Rician_sum=Cor_Rician_sum+real(correlation_rician);
         
end
NorCorrelXYExpOper1=tempIsoXX/MC_max;
NorCorrelXYExpOper2=tempIsoDP/MC_max;  
Cor_Rician=Cor_Rician_sum/MC_max;
 
     

XRange=par.UETX.distance(1,:)/par.wav-par.UETX.distance(1,1)/par.wav;
% XLimRange=[0 5];

% Plot figures
figure
hold on
plot(XRange,NorCorrelXYExpOper1,'-*','Color',[0.89,0.36,0.36],'linewidth',2)  
plot(XRange,Cor_Rician,'-','Color',[0.55,0.51,0.85],'linewidth',2) 
plot(XRange,ClarkeCor,':','Color', [0.36,0.31,0.31],'linewidth',2) 
legend( ['Proposed model ($K=$',num2str(KFactor),')'], ['Rician model($K=$',num2str(KFactor),')'],'Clarke Model','Interpreter','latex')
ylabel('Correlation function','Interpreter','latex')
xlabel('spacing  / $\lambda$','Interpreter','latex')
 
 ylim([-1 1])
 grid on 
xlim([0 5])
 
title(['$K=$', num2str(KFactor), ', $ d_{rt}=$', num2str(abs(dis_rt )),'$\lambda$'],'Interpreter','latex')
 
 

fprintf('end')



