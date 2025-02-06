

clc;
clear;
close all;

%% Basic parameters of cst modeling
Freq=[4,6]; % frequency range(GHz)
freq=(Freq(2)-Freq(1))/2+Freq(1); % operating frequency
wav=3/(0.01*freq); % wavelength
radius=0.01*wav; % radius of diapole
len=0.45*wav; % length of diapole
gap=0.001*wav; % The gap inside the diapole
spacing=0.4*wav;% spacing between adjacent antennas in antenna array
num_x=4; % number of antennas in x-direction
num_y=4; % number of antennas in y-direction
TotalNum=num_x*num_y; % total number of antennas in HMIMO
 
%% CST Initialization
    cst = actxserver('CSTStudio.application.2024'); % load CST (the version is needed for multiple cst verisons)
    mws = invoke(cst, 'NewMWS');% new MWS project 
    app = invoke(mws, 'GetApplicationName');% get the current aap name
    ver = invoke(mws, 'GetApplicationVersion');% get the version num
    invoke(mws, 'FileNew');% new cst file
    path=pwd; %get the path of current .m file 
    filename=['\HMIMO_',num2str(num_x),'_',num2str(num_y),'_',num2str(spacing),'.cst'];
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
invoke(solver,'FrequencyRange',Freq(1),Freq(2));
release(solver)

%% diapole
% component1:solid1
solid=invoke(mws,'Cylinder');
invoke(solid,'Reset');
invoke(solid,'Name','solid1');
invoke(solid,'Component','component1');
invoke(solid,'Material','PEC');
invoke(solid,'OuterRadius',radius);
invoke(solid,'InnerRadius',0);
invoke(solid,'Axis','z');
invoke(solid,'Zrange',-len/2,len/2);
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
invoke(solid,'OuterRadius',radius);
invoke(solid,'InnerRadius',0);
invoke(solid,'Axis','z');
invoke(solid,'Zrange',-gap/2,gap/2);
invoke(solid,'Xcenter',0);
invoke(solid,'Ycenter',0);
invoke(solid,'Create');
release(solid);

% Boolean: Subtract 
solid=invoke(mws,'Solid');
invoke(solid,'Subtract','component1:solid1','component1:solid2')
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
invoke(discreteport,'SetP1','True','-2.3389030995116e-17','0',gap/2);
invoke(discreteport,'SetP2','True','-2.3389030995116e-17','0',gap/2);
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
invoke(transform,'Vector',spacing,0,0); 
invoke(transform,'MultipleObjects','True');
invoke(transform,'Repetitions',num_x-1);
invoke(transform,'Transform','Shape','Translate');
release(transform);
% duplicate the port1
transform=invoke(mws,'Transform');
invoke(transform,'Reset');
invoke(transform,'Name','port1');
invoke(transform,'Vector',spacing,0,0); 
invoke(transform,'MultipleObjects','True');
invoke(transform,'Repetitions',num_x-1);
invoke(transform,'Transform','Port','Translate');
release(transform);
% rename the first component1 (for easiness of subsequent duplication)

for ind_x=1:num_x-1
    solid=invoke(mws,'Solid');
invoke(solid,'Rename',['component1:solid1_',num2str(ind_x)],['component1:solid',num2str(ind_x+1)]);
release(solid); 
end

% Duplicate in y-direction
for ind_x=1:num_x 
% duplicate the compomnent1:solid1
transform=invoke(mws,'Transform');
invoke(transform,'Reset');
invoke(transform,'Name',['component1:solid',num2str(ind_x)]);
invoke(transform,'Vector',0,spacing+len/2,0); 
invoke(transform,'MultipleObjects','True');
invoke(transform,'Repetitions',num_y-1);
invoke(transform,'Transform','Shape','Translate');
release(transform);
% duplicate the port1
transform=invoke(mws,'Transform');
invoke(transform,'Reset');
invoke(transform,'Name',['port',num2str(ind_x)]);
invoke(transform,'Vector',0,spacing+len/2,0); 
invoke(transform,'MultipleObjects','True');
invoke(transform,'Repetitions',num_y-1);
invoke(transform,'Transform','Port','Translate');
release(transform);
end

% rename the component
for ind_x=0:num_x-1
    for ind_y=1:num_y-1
    solid=invoke(mws,'Solid');
    invoke(solid,'Rename',['component1:solid',num2str(ind_x+1),...
         '_',num2str(ind_y)],['component1:solid',num2str(num_x+(num_y-1)*ind_x+ind_y)]);
    release(solid); 
    end
end



%% Substrate material setting
brick=invoke(mws,'Brick');
invoke(brick,'Reset');
invoke(brick,'Name',['solid',num2str(TotalNum+1)]);
invoke(brick,'Component','component1');
invoke(brick,'Material','PEC');
invoke(brick,'Xrange',-1*wav,2*wav);
invoke(brick,'Yrange',-1*wav,2*wav);
invoke(brick,'Zrange',-0.25*wav,-0.25*wav); % thickness is 0
invoke(brick,'Create')
release(brick);

%% Far-field pattern monitor
monitor = invoke(mws, 'Monitor');
farfield_monitor =freq;
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
%% Combination of patterns
solver = invoke(mws, 'Solver');
invoke(solver, 'ResetExcitationModes'); 
invoke(solver, 'SParameterPortExcitation', 'False');
invoke(solver, 'SimultaneousExcitation' ,'True');
invoke(solver, 'SetSimultaneousExcitAutoLabel' ,'True');
invoke(solver,'SetSimultaneousExcitationOffset', 'Phaseshift');
invoke(solver, 'PhaseRefFrequency' ,freq);

% solver = invoke(mws, 'Solver');
% invoke(solver, 'Method','Hexahedral'); 
% invoke(solver, 'StimulationPort ', 'All'); 
% invoke(solver, 'StimulationMode ', 'All'); 

solver = invoke(mws, 'Solver');
invoke(solver, 'Start');
invoke(mws, 'Save');

solver = invoke(mws, 'Solver');
invoke(solver, 'ResetExcitationModes'); 
invoke(solver, 'SParameterPortExcitation', 'False');
invoke(solver, 'SimultaneousExcitation' ,'True');
invoke(solver, 'SetSimultaneousExcitAutoLabel' ,'False');
invoke(solver,'SetSimultaneousExcitationOffset', 'Phaseshift');
invoke(solver, 'PhaseRefFrequency' ,freq);
release(solver);

sCommand = '';
sCommand = [sCommand 'Mesh.SetCreator "High Frequency"  '];

solver = invoke(mws, 'Solver');
invoke(solver, 'Method','Hexahedral');  
invoke(solver, 'CalculationType','TD-S');  
invoke(solver, 'StimulationPort','Selected');  
invoke(solver, 'StimulationMode ', 'All'); 
invoke(solver, 'SteadyStateLimit ', '-40'); 
invoke(solver, 'MeshAdaption ', 'False'); 
release(solver);

solver = invoke(mws, 'Solver');
invoke(solver, 'Start');
invoke(mws, 'Save');

%% Start Simulation
solver = invoke(mws, 'Solver');
invoke(solver, 'Start');
invoke(mws, 'Save');
release(cst);
