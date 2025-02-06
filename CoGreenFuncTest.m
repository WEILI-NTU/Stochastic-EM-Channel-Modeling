function  [CoGreenSP,CoGreenDP,CoGreenTP]=CoGreenFuncTest(par,temp_UETX_dis)
%% Generate the coherent part, i.e., dyadic Green's function
alpha=par.wavenumber*par.UETX.distance;
ChanScal=par.wavenumber*exp(-1i*alpha)./(4*pi*alpha.^3); % N_r*N_s
 
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

CoGreenXX=ChanScal.*PolarFactorXX;
CoGreenYY=ChanScal.*PolarFactorYY;
CoGreenZZ=ChanScal.*PolarFactorZZ;
CoGreenXY=ChanScal.*PolarFactorXY;
CoGreenXZ=ChanScal.*PolarFactorXZ;
CoGreenYZ=ChanScal.*PolarFactorYZ;


CoGreenTP=[CoGreenXX,CoGreenXY,CoGreenXZ;CoGreenXY,CoGreenYY,CoGreenYZ;...
    CoGreenXZ,CoGreenYZ,CoGreenZZ];
CoGreenDP=[CoGreenXX,CoGreenXY,CoGreenXZ*0;CoGreenXY,CoGreenYY,CoGreenYZ*0;...
    CoGreenXZ,CoGreenYZ,CoGreenZZ*0];
CoGreenSP=[CoGreenXX,CoGreenXY*0,CoGreenXZ*0;CoGreenXY,CoGreenYY*0,CoGreenYZ*0;...
    CoGreenXZ,CoGreenYZ*0,CoGreenZZ*0];

end