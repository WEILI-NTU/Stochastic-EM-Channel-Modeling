%==========================================================================
% ==========Stochastic Green's function based channel modeling=============
%==========================================================================
%  Author: Wei Li;
%  Date:   05/07/ 2024,  NTU;
%  Version: V1.0 

%  Note: This function computes the practical and theoretical 
%        capacity region of differnt K-factor

%%========================================================
clc
clear   
%%========================================================  
KFactor=2; % Rician chanel  
users.num=2; 
dis_rt=[200,300]; % distance between the transmitter and receiver
flagUsers='FF-FF'; % Both users are located in "NF-NF" or "FF-FF" or "NF-FF"
%% Simulation setting
MC_max=100; % Monto-Carlo
SNR_range= [ -10  ];
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


users.start.h=0:rec.aperture:users.num*rec.aperture; 
users.start.v=0:rec.aperture:users.num*rec.aperture;  
users.start.z=dis_rt(1:users.num)*par.wav;


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
K_constant=cellfun(@(z)  KFactor/z, RatioLosNlos,'UniformOutput',false); 
flagNLoS=cellfun(@(z) sqrt(1/(1+z)), K_constant,'UniformOutput',false);  
flagLoS=cellfun(@(z) sqrt(1-1/(1+z)), K_constant,'UniformOutput',false); 



varNoise=10.^(-0.1.*SNR_range);
PowerCurrent=eye(tra.totalNum); 


ColorIndex={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],...
    [0.4660 0.6740 0.1880],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330]};


figure
hold on  
for  snrVal=1:length(varNoise)
  
UnCouptempPointA_X={0;0;0;0;0;0};
UnCouptempPointA_Y={0;0;0;0;0;0};
UnCouptempPointC_X={0;0;0;0;0;0};
UnCouptempPointC_Y={0;0;0;0;0;0};

CouptempPointA_X={0;0;0;0;0;0};
CouptempPointA_Y={0;0;0;0;0;0};
CouptempPointC_X={0;0;0;0;0;0};
CouptempPointC_Y={0;0;0;0;0;0};
 

RiciantempPointA_X=0;
RiciantempPointA_Y=0;
RiciantempPointC_X=0;
RiciantempPointC_Y=0;

    for MC=1:MC_max
          NLOSOper= cellfun(@(x,y)  termConstant*termQV*x+ (termConstant*termQV)*sqrt(y).*randn(size(y,1),size(y,2)),IncoGreenExpOper,IncoGreenVarOper ,'UniformOutput',false);
        % K-factor
        NLOSOper= cellfun(@(x,y)  x*y,NLOSOper,flagNLoS ,'UniformOutput',false);
        LOSOper= cellfun(@(x,y)  termConstant*x*y,CoGreenOper,flagLoS ,'UniformOutput',false);
        IsoChOper= cellfun(@(x,y)  x+ y,LOSOper,NLOSOper, 'UniformOutput',false); 

        % Uncoupled channel capacity
        UnCoupleChOper= cellfun(@(z)z,IsoChOper,'UniformOutput',false); 

        UnCoupApointOper_X=cellfun(@(x) sum(log2(1+(eig(...
            x(1,:) * x(1,:)'* pinv(varNoise(snrVal)+x(2,:)* x(2,:)')...
           ) ) ) ),  UnCoupleChOper, 'UniformOutput',false);

        UnCoupApointOper_Y=cellfun(@(x) sum(log2(1+(eig(...
             x(2,:)* x(2,:)'/ varNoise(snrVal))) )),...
            UnCoupleChOper, 'UniformOutput',false);
 
        % Index of point C
        UnCoupCpointOper_X=cellfun(@(x) sum(log2(1+(eig(...
            x(1,:)* x(1,:)'/(varNoise(snrVal)) )) ) ), UnCoupleChOper, 'UniformOutput',false);
        UnCoupCpointOper_Y=cellfun(@(x) sum(log2(1+(eig(...
             x(2,:)* x(2,:)'* pinv(varNoise(snrVal)+x(1,:)*x(1,:)')...
            )) ) ), UnCoupleChOper, 'UniformOutput',false);
         
        % Monto-Carlo Simulations  
        UnCouptempPointA_X=cellfun(@(x,y) x+y, UnCoupApointOper_X,UnCouptempPointA_X, 'UniformOutput',false);
        UnCouptempPointA_Y=cellfun(@(x,y) x+y, UnCoupApointOper_Y,UnCouptempPointA_Y, 'UniformOutput',false);
        UnCouptempPointC_X=cellfun(@(x,y) x+y, UnCoupCpointOper_X,UnCouptempPointC_X, 'UniformOutput',false);
        UnCouptempPointC_Y=cellfun(@(x,y) x+y, UnCoupCpointOper_Y,UnCouptempPointC_Y, 'UniformOutput',false);
 

         %% Rician channel
        LoSChannel= cellfun(@(x)  termConstant*x,CoGreenOper,'UniformOutput',false);
        % Compute the ratio between the LOS and NLOS  
        flagNLoS_Ric=  sqrt(1/(1+KFactor)) ; 
        flagLoS_Ric= sqrt(1-1/(1+KFactor)) ; 
        % To generate a Rayleigh fading channel gain
        h_rayleigh = sqrt(1/2) *(randn(rec.num_x*rec.num_y, ...
            tra.num_x*tra.num_y) + 1i*randn(rec.num_x*rec.num_y,tra.num_x*tra.num_y)); 
        % To generate a Rician fading channel gain
         h_rician = flagLoS_Ric*flagLoS{1} *  (LoSChannel{1}) +flagNLoS_Ric*  h_rayleigh; 
     
         RicianApointOper_X= sum(log2(1+(eig(...
            h_rician(1,:) * h_rician(1,:)'* pinv(varNoise(snrVal)+h_rician(2,:)* h_rician(2,:)')...
           ) ) ) ) ;

        RicianApointOper_Y= sum(log2(1+(eig(...
             h_rician(2,:)* h_rician(2,:)'/ varNoise(snrVal))) )) ;
 
        % Index of point C
        RicianCpointOper_X= sum(log2(1+(eig(...
            h_rician(1,:)* h_rician(1,:)'/(varNoise(snrVal)) )) ) );
        RicianCpointOper_Y= sum(log2(1+(eig(...
             h_rician(2,:)* h_rician(2,:)'* pinv(varNoise(snrVal)+h_rician(1,:)*h_rician(1,:)')...
            )) ) );
         
        % Monto-Carlo Simulations  
        RiciantempPointA_X= RicianApointOper_X+RiciantempPointA_X;
        RiciantempPointA_Y=  RicianApointOper_Y+RiciantempPointA_Y;
        RiciantempPointC_X= RicianCpointOper_X+RiciantempPointC_X;
        RiciantempPointC_Y=  RicianCpointOper_Y+RiciantempPointC_Y;
 
    end
 
    UnCouptempPointA_X=cellfun(@(x)x/MC_max,  UnCouptempPointA_X, 'UniformOutput',false);
     UnCouptempPointA_Y=cellfun(@(x) x/MC_max,  UnCouptempPointA_Y, 'UniformOutput',false);
     UnCouptempPointC_X=cellfun(@(x) x/MC_max, UnCouptempPointC_X, 'UniformOutput',false);
     UnCouptempPointC_Y=cellfun(@(x) x/MC_max, UnCouptempPointC_Y, 'UniformOutput',false);
     
    UnCoupcollect_Apoint_X=UnCouptempPointA_X{1};
    UnCoupcollect_Apoint_Y=UnCouptempPointA_Y{1};
    UnCoupcollect_Cpoint_X=UnCouptempPointC_X{1};
    UnCoupcollect_Cpoint_Y=UnCouptempPointC_Y{1};
    plot([0,UnCoupcollect_Apoint_X,UnCoupcollect_Cpoint_X,UnCoupcollect_Cpoint_X],...
        [UnCoupcollect_Apoint_Y,UnCoupcollect_Apoint_Y,UnCoupcollect_Cpoint_Y,0],...
        '-','color',[0.91,0.47,0.47],'linewidth',2)


    RiciantempPointA_X=  RiciantempPointA_X/MC_max;
    RiciantempPointA_Y=  RiciantempPointA_Y/MC_max;
    RiciantempPointC_X=  RiciantempPointC_X/MC_max;
    RiciantempPointC_Y=   RiciantempPointC_Y/MC_max;
 
    
   

% % plot the parallel and vertical lines of corner points
% plot([0,UnCoupcollect_Cpoint_X],...
%     [UnCoupcollect_Cpoint_Y,UnCoupcollect_Cpoint_Y],...
%     '--*','color',[0.91,0.47,0.47],'linewidth',2)
% plot([UnCoupcollect_Apoint_X,UnCoupcollect_Apoint_X],...
%     [0,UnCoupcollect_Apoint_Y],...
%     '--*','color',[0.91,0.47,0.47],'linewidth',2)

%% Theoretical performance
    [TheoChannel_Ac,TheoChannel_FF,TheoChannel_FF_ap]=MUChApproximate(par,temp_UETX_dis,termConstant,flagLoS{1},flagNLoS{1},users);
 
switch flagUsers
    case 'NF-NF'
        %% Both points in near-field 
        NFApointTheo_X=log2(1+sum(TheoChannel_FF(1,:)) ...
            /(sum(TheoChannel_FF(2,:))+varNoise(snrVal) ) ) ; 
        NFApointTheo_Y=log2(1+sum(TheoChannel_FF(2,:))/varNoise(snrVal) ) ; 
        NFCpointTheo_X=log2(1+sum(TheoChannel_FF(1,:))/varNoise(snrVal) ) ; 
        NFCpointTheo_Y=log2(1+sum(TheoChannel_FF(2,:)) ...
            /(sum(TheoChannel_FF(1,:))+varNoise(snrVal) ) ) ;
        plot([0,NFApointTheo_X,NFCpointTheo_X,NFCpointTheo_X],...
            [NFApointTheo_Y,NFApointTheo_Y,NFCpointTheo_Y,0],...
            '-','color',[0.21,0.69,0.62],'linewidth',2)
        plot([0,RiciantempPointA_X,RiciantempPointC_X,RiciantempPointC_X],...
            [RiciantempPointA_Y,RiciantempPointA_Y,RiciantempPointC_Y,0],...
            '-','color',[0.95,0.55,0.97],'linewidth',2)
        % h(1) = plot(NaN,NaN,'-','Color',[0.91,0.47,0.47],'linewidth',2); 
        % h(2) = plot(NaN,NaN,'-','Color',[0.21,0.69,0.62],'linewidth',2); 
        legend(['Practical (NF-NF, $K$=',num2str(KFactor),')'], ...
            ['Theo. (NF-NF, $K$=',num2str(KFactor),')'],['Rician Channel (NF-NF, $K$=',num2str(KFactor),')'], 'Interpreter','latex')
    case 'FF-FF'
    %% Both points in far-field 
    FFApointTheo_X=log2(1+sum(TheoChannel_FF_ap(1,:)) ...
        /(sum(TheoChannel_FF_ap(2,:))+varNoise(snrVal) ) ) ; 
    FFApointTheo_Y=log2(1+sum(TheoChannel_FF_ap(2,:))/varNoise(snrVal) ) ; 
    FFCpointTheo_X=log2(1+sum(TheoChannel_FF_ap(1,:))/varNoise(snrVal) ) ; 
    FFCpointTheo_Y=log2(1+sum(TheoChannel_FF_ap(2,:)) ...
        /(sum(TheoChannel_FF_ap(1,:))+varNoise(snrVal) ) ) ;
    plot([0,FFApointTheo_X,FFCpointTheo_X,FFCpointTheo_X],...
        [FFApointTheo_Y,FFApointTheo_Y,FFCpointTheo_Y,0],...
        '-','color',[0.45,0.43,0.86],'linewidth',2)
    plot([0,RiciantempPointA_X,RiciantempPointC_X,RiciantempPointC_X],...
            [RiciantempPointA_Y,RiciantempPointA_Y,RiciantempPointC_Y,0],...
            '-','color',[0.95,0.55,0.97],'linewidth',2) 
    legend(['Practical (FF-FF, $K$=',num2str(KFactor),')'], ...
            ['Theo. (FF-FF, $K$=',num2str(KFactor),')'],['Rician Channel (FF-FF,$K$=',num2str(KFactor),')'], 'Interpreter','latex')
    case 'NF-FF'
    %% One near-field user and one far-field user
    NFFFApointTheo_X=log2(1+sum(TheoChannel_FF(1,:)) ...
            /(sum(TheoChannel_FF_ap(2,:))+varNoise(snrVal) ) ) ; 
    NFFFApointTheo_Y=log2(1+sum(TheoChannel_FF_ap(2,:))/varNoise(snrVal) ) ; 
    NFFFCpointTheo_X=log2(1+sum(TheoChannel_FF(1,:))/varNoise(snrVal) ) ; 
    NFFFCpointTheo_Y=log2(1+sum(TheoChannel_FF_ap(2,:)) ...
        /(sum(TheoChannel_FF(1,:))+varNoise(snrVal) ) ) ;
     plot([0,NFFFApointTheo_X,NFFFCpointTheo_X,NFFFCpointTheo_X],...
        [NFFFApointTheo_Y,NFFFApointTheo_Y,NFFFCpointTheo_Y,0],...
        '-','color',[0.29,0.64,0.79],'linewidth',2)
    plot([0,RiciantempPointA_X,RiciantempPointC_X,RiciantempPointC_X],...
            [RiciantempPointA_Y,RiciantempPointA_Y,RiciantempPointC_Y,0],...
            '-','color',[0.95,0.55,0.97],'linewidth',2) 
    legend(['Practical (NF-FF, $K$=',num2str(KFactor),')'], ...
            ['Theo. (NF-FF, $K$=',num2str(KFactor),')'],['Rician Channel (NF-FF,$K$=',num2str(KFactor),')'], 'Interpreter','latex')
    otherwise
        warning('Unexpected users mode')
end



 

 

xlabel('$\mathcal{R}_1$  (bits/s/Hz)','Interpreter','latex');
ylabel('$\mathcal{R}_2$  (bits/s/Hz)','Interpreter','latex');
% legend(Leg_str, 'Interpreter','latex')

title(['$K=$',num2str(KFactor),', Capacity for ',num2str(tra.totalNum),'$\times$', num2str(rec.totalNum),' array',', $d_{rt}=$', num2str(abs(dis_rt)) ],'Interpreter','latex')
grid on
 
end
 
 
fprintf('end')



