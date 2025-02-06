%==========================================================================
% ==========Stochastic Green's function based channel modeling=============
%==========================================================================
%  Author: Wei Li;
%  Date:   04/07/ 2024,  NTU;
%  Version: V1.0 

%  Note: This function computes the autocorrelation function at 
%        the transmitter for the composed channel, LoS, and NLoS 
%        (controlled by flagLOS and flagNLOS) 

% cavity volume and quality factor jointly control the multi-path
% environment

% Please cite the paper properly if you use any part of this code


%%========================================================
clc
clear   
%%======================================================== 
dis_rt=0.4; % distance between the transmitter and receiver
flagLOS=1; % LoS components exist for 1
flagNLOS=1; % NLOS components exist for 1
MC_max=1000; % Monto-Carlo
%%  Set the parameters
par.Freq=[4,6]*1e9; % frequency range(GHz)
par.freq=(par.Freq(2)-par.Freq(1))/2+par.Freq(1); % operating frequency
par.wav=3e8/(par.freq); % wavelength  3/(10*par.freq)?
par.wavenumber=2*pi/par.wav;
alpha= 10;
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
tra.spac=0.05*par.wav;% spacing between adjacent antennas in array
num_Sampling=80;
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

   

% Clarke's spatial autocorrelation function  
ClarkeCor= besselj(0,2*pi*(par.UETX.distance(1,:)-par.UETX.distance(1,1) )/par.wav); 
% ClarkeCor=  sinc((par.UETX.distance(1,:)-par.UETX.distance(1,1) )/par.wav); 
 

 
% Generate isolated channel matrix 
NorCorrelXYExpOper2_pre=[];
NorCorrelXYExpOper3_pre=[];
NorCorrelXYExpOper9_pre=[];


for MC=1:MC_max
    NLOSOper= cellfun(@(x,y)  termConstant*x+ (termConstant)*sqrt(y).*randn(size(y,1),size(y,2)),IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false);
    IsoChOper= cellfun(@(x,y,fl,nl)  termConstant*flagLOS*x+ flagNLOS*y,CoGreenOper,NLOSOper ,'UniformOutput',false);
    IsoChOper= cellfun(@(z) normalize(z,"norm"), IsoChOper,'UniformOutput',false);
    SumIsoChanOper= cellfun(@(z) z'*z, IsoChOper,'UniformOutput',false); 
    SumXYExpOper1= (SumIsoChanOper{2}); % YY-polarized 
    SumXYExpOper2= (SumIsoChanOper{1}+SumIsoChanOper{2} )/2 ; 
    SumXYExpOper3= (SumIsoChanOper{1}+SumIsoChanOper{2}+SumIsoChanOper{3})/3 ;
    SumXYExpOper9= (SumIsoChanOper{1}+SumIsoChanOper{2}+SumIsoChanOper{3}...
        +2*SumIsoChanOper{4}+2*SumIsoChanOper{5}+2*SumIsoChanOper{6})/9 ;
    
    NorCorrelXYExpOper1_pre=   real(SumXYExpOper1(:,1));
    NorCorrelXYExpOper2_pre=  real(SumXYExpOper2(1,:));
    NorCorrelXYExpOper3_pre=  real(SumXYExpOper3(1,:));
    NorCorrelXYExpOper9_pre=  real(SumXYExpOper9(1,:));

% NorCorrelXYExpOper1_pre=  real(SumXYExpOper1(:,1)/max(real(SumXYExpOper1(:,1))));
%     NorCorrelXYExpOper2_pre=  real(SumXYExpOper2(1,:)/max(real(SumXYExpOper2(1,:))));
%     NorCorrelXYExpOper3_pre=  real(SumXYExpOper3(1,:)/max(real(SumXYExpOper3(1,:))));
%     NorCorrelXYExpOper9_pre=  real(SumXYExpOper9(1,:)/max(real(SumXYExpOper9(1,:))));

    if MC==1
        NorCorrelXYExpOper1=NorCorrelXYExpOper1_pre;
        NorCorrelXYExpOper2=NorCorrelXYExpOper2_pre;
        NorCorrelXYExpOper3=NorCorrelXYExpOper3_pre;
        NorCorrelXYExpOper9=NorCorrelXYExpOper9_pre;
    else
        NorCorrelXYExpOper1=NorCorrelXYExpOper1+NorCorrelXYExpOper1_pre;
        NorCorrelXYExpOper2=NorCorrelXYExpOper2+NorCorrelXYExpOper2_pre;
        NorCorrelXYExpOper3=NorCorrelXYExpOper3+NorCorrelXYExpOper3_pre;
        NorCorrelXYExpOper9=NorCorrelXYExpOper9+NorCorrelXYExpOper9_pre;
    end
end
NorCorrelXYExpOper1=NorCorrelXYExpOper1/MC_max;
NorCorrelXYExpOper2=NorCorrelXYExpOper2/MC_max;
NorCorrelXYExpOper3=NorCorrelXYExpOper3/MC_max;
NorCorrelXYExpOper9=NorCorrelXYExpOper9/MC_max;

NorCorrelXYExpOper1=  NorCorrelXYExpOper1/(max(NorCorrelXYExpOper1));
    NorCorrelXYExpOper2=  NorCorrelXYExpOper2/(max(NorCorrelXYExpOper2));
    NorCorrelXYExpOper3= NorCorrelXYExpOper3/(max(NorCorrelXYExpOper3));
    NorCorrelXYExpOper9=  NorCorrelXYExpOper9/(max(NorCorrelXYExpOper9));
     

XRange=par.UETX.distance(1,:)/par.wav-par.UETX.distance(1,1)/par.wav;
XLimRange=[0 5];

% Plot figures for all polarization states
% figure
% hold on
% plot(XRange,NorCorrelXYExpOper1,'sr','linewidth',2) 
% plot(XRange,NorCorrelXYExpOper2,'xb','linewidth',2) 
% plot(XRange,NorCorrelXYExpOper3,'o-m','linewidth',2) 
% plot(XRange,NorCorrelXYExpOper9,'+','linewidth',2) 
% plot(XRange,ClarkeCor,'-g','Color',[0.25 0.80 0.54],'linewidth',2) 
% legend( 'YY','$(H_{xx}^2+H_{yy}^2)/2$','$(H_{xx}^2+H_{yy}^2+H_{zz}^2)/3$','co/cross-polarization channels','Clarke Model','Interpreter','latex')
% ylabel('Correlation function','Interpreter','latex')
% xlabel('spacing  / $\lambda$','Interpreter','latex')
%  % xlim(XLimRange)
%  ylim([-1 1])
%  grid on 

 % Plot figures for single polarization state 
figure
hold on
plot(XRange,NorCorrelXYExpOper1,'-','Color',[0.86,0.27,0.27],'linewidth',2)  
plot(XRange,ClarkeCor,':o','Color', [0.36,0.31,0.31],'linewidth',2) 
legend( 'Proposed model','Clarke Model','Interpreter','latex')
ylabel('Correlation function','Interpreter','latex')
xlabel('spacing  / $\lambda$','Interpreter','latex')
 % xlim(XLimRange)
 ylim([-1 1])
 grid on 
xlim(XLimRange)
 
title(['fast,', num2str([flagLOS,flagNLOS]), ', $ d_{rt}=$', num2str(abs(dis_rt )),'$\lambda$, Q is ', num2str(par.Qfactor),', V is ', num2str(par.cavity/par.wav),'$\lambda$'],'Interpreter','latex')
 
 

fprintf('end')



