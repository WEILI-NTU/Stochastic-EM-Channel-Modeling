function  [ChanPol]=DP_ChannelGen(rec,par,temp_UETX_dis)
% This function generates the near field channel
 
%% Generate the channel using Green's function
alpha=2*pi/par.wav*par.UETX.distance;
ChanScalFirstTerm=1j*par.wavenumber*120*pi*exp(-1i*alpha)./(4*pi*alpha.^3); % N_r*N_s

tempSincX=par.wavenumber*temp_UETX_dis.X*rec.spac...
         ./(2*par.UETX.distance);
ChanScalSecondTerm=sin(tempSincX)./tempSincX;
tempSincY=par.wavenumber*temp_UETX_dis.Y*rec.spac...
    ./(2*par.UETX.distance);
ChanScalThirdTerm=sin(tempSincY)./tempSincY;
sincTerm=ChanScalSecondTerm.*ChanScalThirdTerm;

sincTerm(isnan(sincTerm))=1;

ChanScal=ChanScalFirstTerm.*sincTerm; 
 


% The direction of each tra-rec pair
temp_direc_x= -temp_UETX_dis.X./par.UETX.distance;
temp_direc_y= -temp_UETX_dis.Y./par.UETX.distance;
temp_direc_z= temp_UETX_dis.Z./par.UETX.distance;

direc.XX=temp_direc_x.*temp_direc_x;
direc.YY=temp_direc_y.*temp_direc_y;
direc.ZZ=temp_direc_z.*temp_direc_z;
direc.XY=temp_direc_x.*temp_direc_y;
direc.XZ=temp_direc_x.*temp_direc_z;
direc.YZ=temp_direc_y.*temp_direc_z;


ChanPolarFactorFirst=alpha.^2-1-1i.*alpha;
ChanPolarFactorSecond=3-alpha.^2+3*1i.*alpha;

PolarFactorXX=ChanPolarFactorFirst+ChanPolarFactorSecond.*direc.XX;
PolarFactorYY=ChanPolarFactorFirst+ChanPolarFactorSecond.*direc.YY;
PolarFactorZZ=ChanPolarFactorFirst+ChanPolarFactorSecond.*direc.ZZ;
PolarFactorXY= ChanPolarFactorSecond.*direc.XY;
PolarFactorXZ= ChanPolarFactorSecond.*direc.XZ;
PolarFactorYZ= ChanPolarFactorSecond.*direc.YZ;

ChanPol.XX=ChanScal.*PolarFactorXX;
ChanPol.XY=ChanScal.*PolarFactorXY;
ChanPol.XZ=ChanScal.*PolarFactorXZ;
ChanPol.YY=ChanScal.*PolarFactorYY;
ChanPol.YZ=ChanScal.*PolarFactorYZ;
ChanPol.ZZ=ChanScal.*PolarFactorZZ;

end