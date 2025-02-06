function CouplingMatrixGen(par,tra)
% This function generates the impedance matrix and compute the mutual
% coupling matrix
 

%% Basic parameters of cst modeling
% par.Freq=[4,6]; % frequency range(GHz)
% % freq=(par.Freq(2)-par.Freq(1))/2+par.Freq(1); % operating frequency
% % par.wav=3/(0.01*par.freq); % wavelength
% par.rad=0.01*par.wav; % radius of diapole
% par.len=0.43*par.wav; % length of diapole
% par.gap=0.01*par.wav; % The gap inside the diapole
% tra.spac=0.4*par.wav;% spacing between adjacent antennas in antenna array
% tra.num_x=2; % number of antennas in x-direction
% tra.num_y=2; % number of antennas in y-direction
% tra.totalNum=tra.num_x*tra.num_y; % total number of antennas in HMIMO
 
%% CST Initialization
    cst = actxserver('CSTStudio.application.2024'); % load CST (the version is needed for multiple cst verisons)
    mws = invoke(cst, 'NewMWS');% new MWS project 
    app = invoke(mws, 'GetApplicationName');% get the current aap name
    ver = invoke(mws, 'GetApplicationVersion');% get the version num
    invoke(mws, 'FileNew');% new cst file
    path=pwd; %get the path of current .m file 
    filename=['\HMIMO_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav*10),'_',num2str(par.freq),'.cst'];
    fullname=[path filename];
    invoke(mws, 'SaveAs', fullname, 'True');%True denote the result saved currently
    invoke(mws, 'DeleteResults'); % delete the previous result

%% Global Unit Initialization 
units = invoke(mws, 'Units');
invoke(units, 'Geometry', 'mm');
invoke(units, 'Frequency', 'ghz');
invoke(units, 'Time', 'ns');
invoke(units, 'TemperatureUnit', 'kelvin');

%% Frequency setting
solver=invoke(mws,'Solver');
invoke(solver,'FrequencyRange',par.Freq(1),par.Freq(2));
release(solver)

%% diapole
% component1:solid1
solid=invoke(mws,'Cylinder');
invoke(solid,'Reset');
invoke(solid,'Name','solid1');
invoke(solid,'Component','component1');
invoke(solid,'Material','PEC');
invoke(solid,'OuterRadius',par.rad);
invoke(solid,'InnerRadius',0);
invoke(solid,'Axis','z');
invoke(solid,'Zrange',-par.len/2,par.len/2);
invoke(solid,'Xcenter',0);
invoke(solid,'Ycenter',0);
invoke(solid,'Create');
release(solid);

% component1:solid2
solid=invoke(mws,'Cylinder');
invoke(solid,'Reset');
invoke(solid,'Name','solid2');
invoke(solid,'Component','component1');
invoke(solid,'Material','PEC');
invoke(solid,'OuterRadius',par.rad);
invoke(solid,'InnerRadius',0);
invoke(solid,'Axis','z');
invoke(solid,'Zrange',-par.gap/2,par.gap/2);
invoke(solid,'Xcenter',0);
invoke(solid,'Ycenter',0);
invoke(solid,'Create');
release(solid);

% Boolean: Subtract 
solid=invoke(mws,'Solid');
invoke(solid,'Subtract','component1:solid1','component1:solid2');
release(solid);

% Pick center point at one side
pick=invoke(mws,'Pick');
invoke(pick,'PickCenterpointFromId','component1:solid1','1');
release(pick);

% Pick center point at the other side
pick=invoke(mws,'Pick');
invoke(pick,'PickCenterpointFromId','component1:solid1','3');
release(pick);

% Define discrete port
discreteport=invoke(mws,'DiscretePort');
invoke(discreteport,'Reset');
invoke(discreteport,'PortNumber','1');
invoke(discreteport,'Impedance','50.0');
invoke(discreteport,'Monitor','True');
invoke(discreteport,'Wire','');
invoke(discreteport,'Position','end1'); 
invoke(discreteport,'Radius','0.0');
invoke(discreteport,'SetP1','True','-2.3389030995116e-17','0',par.gap/2);
invoke(discreteport,'SetP2','True','-2.3389030995116e-17','0',par.gap/2);
invoke(discreteport,'InvertDirection','False');
invoke(discreteport,'LocalCoordinates','False'); 
invoke(discreteport,'Create');
release(discreteport);

%% Rotate the diapole
% rotate the compomnent1:solid1
transform=invoke(mws,'Transform');
invoke(transform,'Name','component1:solid1');
invoke(transform,'Origin','Free');
invoke(transform,'Center',0,0,0);
invoke(transform,'Angle','90','0','0');
invoke(transform,'MultipleObjects','False');
invoke(transform,'Repetitions','1');
invoke(transform,'Transform','Shape','Rotate');
release(transform);
% rotate the port1
transform=invoke(mws,'Transform');
invoke(transform,'Name','port1');
invoke(transform,'Origin','Free');
invoke(transform,'Center',0,0,0);
invoke(transform,'Angle','90','0','0');
invoke(transform,'MultipleObjects','False');
invoke(transform,'Repetitions','1');
invoke(transform,'Transform','Port','Rotate');
release(transform);

%% Duplicate diapoles 
% Duplicate in x-direction
% duplicate the compomnent1:solid1
transform=invoke(mws,'Transform');
invoke(transform,'Reset');
invoke(transform,'Name','component1:solid1');
invoke(transform,'Vector',tra.spac+par.rad,0,0); 
invoke(transform,'MultipleObjects','True');
invoke(transform,'Repetitions',tra.num_x-1);
invoke(transform,'Transform','Shape','Translate');
release(transform);
% duplicate the port1
transform=invoke(mws,'Transform');
invoke(transform,'Reset');
invoke(transform,'Name','port1');
invoke(transform,'Vector',tra.spac+par.rad,0,0); 
invoke(transform,'MultipleObjects','True');
invoke(transform,'Repetitions',tra.num_x-1);
invoke(transform,'Transform','Port','Translate');
release(transform);
% rename the first component1 (for easiness of subsequent duplication)

for ind_x=1:tra.num_x-1
    solid=invoke(mws,'Solid');
invoke(solid,'Rename',['component1:solid1_',num2str(ind_x)],['component1:solid',num2str(ind_x+1)]);
release(solid); 
end

% Duplicate in y-direction
for ind_x=1:tra.num_x 
% duplicate the compomnent1:solid1
transform=invoke(mws,'Transform');
invoke(transform,'Reset');
invoke(transform,'Name',['component1:solid',num2str(ind_x)]);
invoke(transform,'Vector',0,tra.spac+par.len/2,0); 
invoke(transform,'MultipleObjects','True');
invoke(transform,'Repetitions',tra.num_y-1);
invoke(transform,'Transform','Shape','Translate');
release(transform);
% duplicate the port1
transform=invoke(mws,'Transform');
invoke(transform,'Reset');
invoke(transform,'Name',['port',num2str(ind_x)]);
invoke(transform,'Vector',0,tra.spac+par.len/2,0); 
invoke(transform,'MultipleObjects','True');
invoke(transform,'Repetitions',tra.num_y-1);
invoke(transform,'Transform','Port','Translate');
release(transform);
end

% rename the component
for ind_x=0:tra.num_x-1
    for ind_y=1:tra.num_y-1
    solid=invoke(mws,'Solid');
    invoke(solid,'Rename',['component1:solid',num2str(ind_x+1),...
         '_',num2str(ind_y)],['component1:solid',num2str(tra.num_x+(tra.num_y-1)*ind_x+ind_y)]);
    release(solid); 
    end
end



%% Substrate material setting
brick=invoke(mws,'Brick');
invoke(brick,'Reset');
invoke(brick,'Name',['solid',num2str(tra.totalNum+1)]);
invoke(brick,'Component','component1');
invoke(brick,'Material','PEC');
invoke(brick,'Xrange',-1*par.wav,3*par.wav);
invoke(brick,'Yrange',-1*par.wav,3*par.wav);
invoke(brick,'Zrange',-0.25*par.wav,-0.25*par.wav); % thickness is 0
invoke(brick,'Create')
release(brick);

%% Far-field pattern monitor
monitor = invoke(mws, 'Monitor');
farfield_monitor =par.freq;
for i = 1:length(farfield_monitor)
    Str_name = ['Farfield (f=',num2str(farfield_monitor(i)),')'];
    invoke(monitor, 'Reset');
    invoke(monitor, 'Name', Str_name);
    invoke(monitor, 'Dimension', 'Volume');
    invoke(monitor, 'Domain', 'Frequency');
    invoke(monitor, 'FieldType', 'Farfield');
    invoke(monitor, 'Frequency', farfield_monitor(i));
    invoke(monitor, 'Create');
end
%% Farfield plot options
farfieldplot = invoke(mws,'FarfieldPlot'); 
invoke(farfieldplot, 'Plottype','3D');
release(farfieldplot);
%% Post-Processing
PostProcess1D = invoke(mws,'PostProcess1D');
invoke(PostProcess1D,'ActivateOperation','yz-matrices','true');
%% Excite separately
solver = invoke(mws, 'Solver');
invoke(solver, 'Start');
invoke(mws, 'Save');

%% Excite simultaneously
solver = invoke(mws, 'Solver');
invoke(solver, 'StimulationPort','Selected'); 
invoke(solver, 'ResetExcitationModes'); 
invoke(solver, 'SParameterPortExcitation', 'False');
invoke(solver, 'SimultaneousExcitation' ,'True');
invoke(solver, 'SetSimultaneousExcitAutoLabel' ,'False');
invoke(solver,'SetSimultaneousExcitationLabel', '1[1.0,0.0]+2[1.0,0.0]+3[1.0,0.0]+4[1.0,0.0],[5]');
invoke(solver,'SetSimultaneousExcitationOffset', 'Phaseshift');
invoke(solver,'ExcitationPortMode', '1', '1', '1.0', '0.0', 'default', 'True' );
invoke(solver,'ExcitationPortMode', '2', '1', '1.0', '0.0', 'default', 'True' );
invoke(solver,'ExcitationPortMode', '3', '1', '1.0', '0.0', 'default', 'True' );
invoke(solver,'ExcitationPortMode', '4', '1', '1.0', '0.0', 'default', 'True' );
invoke(solver, 'PhaseRefFrequency' ,par.freq);
invoke(solver, 'Start');
invoke(mws, 'Save');
release(solver);
 


% %% Export S parameters 
% fid=fopen(['SParameters_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav),'_',num2str(par.freq),'.txt'],'w+');
% fclose(fid);
% exportpath=which(['Zmatrix_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav),'_',num2str(par.freq),'.txt']);
% SelectTreeItem = invoke(mws,'SelectTreeItem','1D Results\S-Parameters');
% ASCIIExport = invoke(mws,'ASCIIExport');
% invoke(ASCIIExport,'Reset');
% invoke(ASCIIExport,'SetVersion','2024');
% invoke(ASCIIExport,'FileName',exportpath);
% invoke(ASCIIExport,'Execute');

%% Export Z parameters 
fid=fopen(['ZMatrix_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav),'_',num2str(par.freq),'.txt'],'w+');
fclose(fid);
exportpath_Z=which(['Zmatrix_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav),'_',num2str(par.freq),'.txt']);
invoke(mws,'SelectTreeItem','1D Results\Z Matrix');
ASCIIExport = invoke(mws,'ASCIIExport');
invoke(ASCIIExport,'Reset'); 
invoke(ASCIIExport,'FileName',exportpath_Z);
invoke(ASCIIExport,'Execute');
release(ASCIIExport) 
 


%% Far-field of the array
% fid=fopen(['Farfield_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav),'_',num2str(par.freq),'.txt'],'w+');
% fclose(fid);
% exportpath_FarField=which(['Farfield_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav),'_',num2str(par.freq),'.txt']);
% invoke(mws,'SelectTreeItem',['Farfields\Farfield (f=',num2str(par.freq),') [1[1.0,0.0]+2[1.0,0.0]+3[1.0,0.0]+4[1.0,0.0],[5]]']);
% ASCIIExport =invoke(mws,'ASCIIExport');
% invoke(ASCIIExport,'Reset'); 
% invoke(ASCIIExport,'FileName',exportpath_FarField);
% invoke(ASCIIExport,'Execute');
% 
% fid=fopen(['Farfield_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav),'_',num2str(par.freq),'.txt'],'w+');
% fclose(fid);
% exportpath_FarField=which(['Farfield_',num2str(tra.num_x),'_',num2str(tra.num_y),'_',num2str(tra.spac/par.wav),'_',num2str(par.freq),'.txt']);
%  port=1;
% SelectTreeItem = invoke(mws,'SelectTreeItem',char(strcat('Farfields\Farfield',num2str(par.freq),strcat({' '}, strcat('[',num2str(port),']'))))); 
% invoke(mws,'SelectTreeItem','Farfields\Farfield (f=5) [1]');
% ASCIIExport =invoke(mws,'ASCIIExport');
% invoke(ASCIIExport,'Reset'); 
% invoke(ASCIIExport,'FileName',exportpath_FarField);
% invoke(ASCIIExport,'Execute');
% 
% SelectTreeItem = invoke(mws,'SelectTreeItem',char(strcat('Farfields\Farfield',num2str(par.freq),strcat({' '}, ' [1[1.0,0.0]+2[1.0,0.0]+3[1.0,0.0]+4[1.0,0.0],[5]]'))));



end