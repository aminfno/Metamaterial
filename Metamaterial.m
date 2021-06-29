
function varargout = Metamaterial(varargin)
%% METAMATERIAL MATLAB code for Metamaterial.fig
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Metamaterial_OpeningFcn, ...
    'gui_OutputFcn',  @Metamaterial_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function Metamaterial_OpeningFcn(hObject, ~, handles, varargin)
load data.mat counter DesignPar Newpar

if counter == 0
else
    DesignPar=Newpar;
    
   %DesignPar=[E,v,H,T,t1,t2,S1,S2,S3,S4,S5,S6,theta1,theta2,strain,nn];
      counter =0;
   save data counter -append  
end
     set(handles.E,'String',DesignPar(1)); 
     set(handles.Poisson,'String',DesignPar(2)); 
     set(handles.H,'String',DesignPar(3)); 
     set(handles.T,'String',DesignPar(4)); 
     set(handles.t1,'String',DesignPar(5)); 
     set(handles.t2,'String',DesignPar(6)); 
     set(handles.S1,'String',DesignPar(7)); 
     set(handles.S2,'String',DesignPar(8)); 
     set(handles.S3,'String',DesignPar(9)); 
     set(handles.S4,'String',DesignPar(10)); 
     set(handles.S5,'String',DesignPar(11)); 
     set(handles.S6,'String',DesignPar(12)); 
     set(handles.theta1,'String',DesignPar(13)); 
     set(handles.theta2,'String',DesignPar(14)); 
     set(handles.strain,'String',DesignPar(15)); 
     set(handles.nn,'String',DesignPar(16)); 
    
    
load data.mat  logo1 reference labels1 labels2 Young  %Load images data


%% GUI images
opengl hardware
axes(handles.LOGO) % Load the Flexible Research Group logo and show it on the GUI
imshow(logo1)
grid off

set(handles.LOGO, 'Box', 'off')
handles.LOGO.XRuler.Axle.LineStyle = 'none';
handles.LOGO.YRuler.Axle.LineStyle = 'none';

axes(handles.Young) % Load the Young modulus and poissons ratio images
imshow(Young)
grid off

set(handles.Young, 'Box', 'off')
handles.Young.XRuler.Axle.LineStyle = 'none';
handles.Young.YRuler.Axle.LineStyle = 'none';

axes(handles.reference) % Load the reference image
imshow(reference);
grid off

set(handles.reference, 'Box', 'off')
handles.reference.XRuler.Axle.LineStyle = 'none';
handles.reference.YRuler.Axle.LineStyle = 'none';
axes(handles.labels1) % Load the t_1 image
imshow(labels1)
grid off

set(handles.labels1, 'Box', 'off')
handles.labels1.XRuler.Axle.LineStyle = 'none';
handles.labels1.YRuler.Axle.LineStyle = 'none';

axes(handles.labels2) % Load the lable2 image
imshow(labels2)
grid off

set(handles.labels2, 'Box', 'off')
handles.labels2.XRuler.Axle.LineStyle = 'none';
handles.labels2.YRuler.Axle.LineStyle = 'none';

handles.output = hObject;
% Update handles structure
guidata(hObject, handles);


function varargout = Metamaterial_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;



function E_Callback(hObject, eventdata, handles)
function E_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Poisson_Callback(hObject, eventdata, handles)
function Poisson_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function H_Callback(hObject, eventdata, handles)
function H_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function T_Callback(hObject, eventdata, handles)

function T_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)

E=str2double(get(handles.E, 'string'));
v=str2double(get(handles.Poisson, 'string'));
H=str2double(get(handles.H, 'string'));
T=str2double(get(handles.T, 'string'));
t1=str2double(get(handles.t1, 'string'));
t2=str2double(get(handles.t2, 'string'));
S1=str2double(get(handles.S1, 'string'));
S2=str2double(get(handles.S2, 'string'));
S3=str2double(get(handles.S3, 'string'));
S4=str2double(get(handles.S4, 'string'));
S5=str2double(get(handles.S5, 'string'));
S6=str2double(get(handles.S6, 'string'));
theta1=str2double(get(handles.theta1, 'string'));
theta2=str2double(get(handles.theta2, 'string'));
strain=str2double(get(handles.strain, 'string'));
nn=str2double(get(handles.nn, 'string'));
DesignPar=[E,v,H,T,t1,t2,S1,S2,S3,S4,S5,S6,theta1,theta2,strain,nn];
Newpar=DesignPar;
save data Newpar -append

%% Generating error string
if sum(isnan(DesignPar))>0
    
    str=sprintf('Please check the following field(s): \n [');
    if isnan(DesignPar(1))==1;str=sprintf('%s E,',str);end
    if isnan(DesignPar(2))==1;str=sprintf('%s v,',str);end
    if isnan(DesignPar(3))==1;str=sprintf('%s H,',str);end
    if isnan(DesignPar(4))==1;str=sprintf('%s T,',str);end
    if isnan(DesignPar(5))==1;str=sprintf('%s t1,',str);end
    if isnan(DesignPar(6))==1;str=sprintf('%s t2,',str);end
    if isnan(DesignPar(7))==1;str=sprintf('%s S1,',str);end
    if isnan(DesignPar(8))==1;str=sprintf('%s S2,',str);end
    if isnan(DesignPar(9))==1;str=sprintf('%s S3,',str);end
    if isnan(DesignPar(10))==1;str=sprintf('%s S4,',str);end
     if isnan(DesignPar(11))==1;str=sprintf('%s S5,',str);end
    if isnan(DesignPar(12))==1;str=sprintf('%s S6,',str);end
    if isnan(DesignPar(13))==1;str=sprintf('%s theta1,',str);end
    if isnan(DesignPar(14))==1;str=sprintf('%s theta2,',str);end
    if isnan(DesignPar(15))==1;str=sprintf('%s Strain amplitude,',str);end
    if isnan(DesignPar(16))==1;str=sprintf('%s Accuracy,',str);end
    
    % Display the final string in the msgbox
    str=strcat(str(1:end-1),' ]');
    errordlg(str,'Error')
    
    return
end


%msg= msgbox('.           Simulation is running ...            .','Please Wait');
global sfig
sfig = waitbar(0,'Please wait...','Name','Please Wait','WindowStyle','modal');
status(1);
tic
Output=DesignTool(DesignPar);%%status(2) updated in DesignTool.m
toc
status(3);%%Update Status
images=cell(size(Output)-1);
for n= 1:size(Output,2)-1
    % Capture the plot as an image
   im =Output{n+1}.cdata(:,:,1);
    images{n}= im;
    if n==150;status(4);end
end
save data images -append
%close(msg);
status(5);%%close status
close(Metamaterial);
Metamaterial_results




function strain_Callback(hObject, eventdata, handles)

function strain_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t1_Callback(hObject, eventdata, handles)

function t1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function t2_Callback(hObject, eventdata, handles)

function t2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function S1_Callback(hObject, eventdata, handles)

function S1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function S2_Callback(hObject, eventdata, handles)

function S2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function S3_Callback(hObject, eventdata, handles)

function S3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function S4_Callback(hObject, eventdata, handles)

function S4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta1_Callback(hObject, eventdata, handles)

function theta1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function theta2_Callback(hObject, eventdata, handles)

function theta2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function E_ButtonDownFcn(hObject, eventdata, handles)
function E_KeyPressFcn(hObject, eventdata, handles)

function status(svar)
global sfig
if svar==1
waitbar(.05,sfig,'Creating the design...');
elseif svar==2
waitbar(.12,sfig,'Processing the design...');
elseif svar==3
waitbar(.4,sfig,'Creating output animation...');
elseif svar==4
    waitbar(.8,sfig,'Finalizing...');
else
close(sfig)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis Code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Output=DesignTool(DesignPar)
% This function must be run to launch the design tool
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Values to input in the GUI
% Set material properties
% E=0.41e9; %Young's modulus of Teflon
E=DesignPar(1);
% v=0.46; %Poisson's ratio of Teflon
v=DesignPar(2);
%  
% % Set geometric parameters for the lattice's cell as defined in the figure
% H=100e-3; %units in meters
H=DesignPar(3);
% T=10e-3; %units in meters
T=DesignPar(4);
% t1=7e-3; %units in meters
t1=DesignPar(5);
% t2=3e-3; %units in meters
t2=DesignPar(6);
% S1=15e-3; %units in meters
S1=DesignPar(7);
% S2=20e-3; %units in meters
S2=DesignPar(8);
% S3=5e-3; %units in meters
S3=DesignPar(9);
% S4=7e-3; %units in meters
S4=DesignPar(10);
S5= DesignPar(11); 
S6= DesignPar(12); 
% theta1=27*pi/180; %units in radians
theta1=DesignPar(13);
% theta2=27*pi/180; %units in radians
theta2=DesignPar(14);
%  
% Specify Loading parameter
% Strain amplitude of oscillating load
strainamplitude=DesignPar(15); %absolute value of the amplitude strain load is 20%

nn=int16(DesignPar(16)/4)*4; % nn points over the entire cycle

% the following should be defined as a function of everthing before
L5=(2*sqrt((((H-(3*T))/(2*cos(theta2)))^2)-(((((H-(3*T))/2)*tan(theta2))+S4)^2)))+(3*T)+(2*t1/cos(theta1));
L6=(3*t1)+T+((H-(3*t1)-T)/cos(theta1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are the dimensions of the parameterized cell shown in the original image
% H=100e-3; %units in meters
% T=9e-3; %units in meters
% t1=3e-3; %units in meters
% t2=6e-3; %units in meters
% S1=20e-3; %units in meters
% S2=16e-3; %units in meters
% S3=10e-3; %units in meters
% S4=20e-3; %units in meters
% theta1=20*pi/180; %units in radians
% theta2=20*pi/180; %units in radians

% Internally set these parameters
% Time to load a single cycle
hr=1*60*60; % 1 hour in units of seconds
% How many points do they want to plot within that cycle

% Specify accuracy term
HH=5; %The larger this is the greater the accuracy but the slower the computation.
% This value doesn't affect the results so we set it equal to T since that is what it will always be when the design is 3D.
W=T; %out-of-plane thickness
% Convert the strain amplitude from a percentage to a ratio
la=strainamplitude/100; 

% Calculate Shear modulus
G=E/(2*(1+v)); %Shear modulus

% Determine how large display window should be. This is internally set.
Horiz=(H/2)+S4+S5; %nnew
Vertic=(H/2)+(((H/2)-(1.5*t1)-(T/2))*tan(theta1))+(t1/cos(theta1))+S6; %nnew
if(Horiz >= Vertic)
    dim=(1+la)*Horiz;
else
    dim=(1+la)*Vertic;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=nn/4; % n points over each quarter of the cycle
tinc=hr/nn;
time=double(0:tinc:hr);
strainload=la*sin((2*pi/hr)*time);
strainPoisson=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the two end points of each of the 20 elements in their initial
% location according to the labeled element and body numbers in the figure.
% First number corresponds to the element number, second number corresponds
% to the body number that the element attaches to. They are separated by an
% s

p1s2=[-H/2, -H/2, 0]+[(1.5*t1), 0, 0]+[-((t1/2)*sin(theta1)), -((t1/2)*cos(theta1)), 0];
p1s0=[-T/2, -(H/2)-(((H/2)-(1.5*t1)-(T/2))*tan(theta1)), 0]+[0, -((t1/2)/cos(theta1)), 0];

p2s1=[-p1s2(1,1), p1s2(1,2), 0];
p2s0=[-p1s0(1,1), p1s0(1,2), 0];

p3s2=[-((H/2)-T-(t2/cos(theta2))), -(H/2)+(T/2), 0];
p3s16=[p3s2(1,1)+S1, p3s2(1,2), 0];

p4s1=[-p3s2(1,1), p3s2(1,2), 0];
p4s15=[-p3s16(1,1), p3s16(1,2), 0];

p5s2=[-(H/2)+(T/2), -(H/2)+T, 0];
p5s18=[p5s2(1,1), p5s2(1,2)+S2, 0];

p6s2=[-((H/2)-T-((t2/2)/cos(theta2))), -(H/2)+T, 0];
p6s6=[-(H/2)+T+(((H/2)-T-(T/2))*tan(theta2))+((t2/2)/cos(theta2)), -T/2, 0];

p7s1=[-p6s2(1,1), p6s2(1,2), 0];
p7s5=[-p6s6(1,1), p6s6(1,2), 0];

p8s1=[-p5s2(1,1), p5s2(1,2), 0];
p8s17=[-p5s18(1,1), p5s18(1,2), 0];

p9s14=[-(H/2)-S4, 0, 0];
p9s6=[-(H/2)+T+(((H/2)-T-(T/2))*tan(theta2)), 0, 0];

p10s6=[-(H/2)+T+(((H/2)-T-(T/2))*tan(theta2))+(t2/cos(theta2)), 0, 0];
p10s12=[p10s6(1,1)+S3, 0, 0];

p11s5=[-p10s6(1,1), 0, 0];
p11s11=[-p10s12(1,1), 0, 0];

p12s13=[-p9s14(1,1), 0, 0];
p12s5=[-p9s6(1,1), 0, 0];

p13s4=[p6s2(1,1), -p6s2(1,2), 0];
p13s6=[p6s6(1,1), -p6s6(1,2), 0];

p14s3=[p7s1(1,1), -p7s1(1,2), 0];
p14s5=[p7s5(1,1), -p7s5(1,2), 0];

p15s4=[p5s2(1,1), -p5s2(1,2), 0];
p15s10=[p5s18(1,1), -p5s18(1,2), 0];

p16s3=[p8s1(1,1), -p8s1(1,2), 0];
p16s9=[p8s17(1,1), -p8s17(1,2), 0];

p17s4=[p3s2(1,1), -p3s2(1,2), 0];
p17s8=[p3s16(1,1), -p3s16(1,2), 0];

p18s3=[p4s1(1,1), -p4s1(1,2), 0];
p18s7=[p4s15(1,1), -p4s15(1,2), 0];

p19s4=[p1s2(1,1), -p1s2(1,2), 0];
p19s19=[p1s0(1,1), -p1s0(1,2), 0];

p20s3=[p2s1(1,1), -p2s1(1,2), 0];
p20s19=[p2s0(1,1), -p2s0(1,2), 0];

%This matrix is added to in the tension scenario
Element=[p1s2, p1s0, t1, 2, 0;
    p2s1, p2s0, t1, 1, 0;
    p3s2, p3s16, T, 2, 16;
    p4s1, p4s15, T, 1, 15;
    p5s2, p5s18, T, 2, 18;
    p6s2, p6s6, t2, 2, 6;
    p7s1, p7s5, t2, 1, 5;
    p8s1, p8s17, T, 1, 17;
    p9s14, p9s6, T, 14, 6;
    p10s6, p10s12, T, 6, 12;
    p11s5, p11s11, T, 5, 11;
    p12s13, p12s5, T, 13, 5;
    p13s4, p13s6, t2, 4, 6;
    p14s3, p14s5, t2, 3, 5;
    p15s4, p15s10, T, 4, 10;
    p16s3, p16s9, T, 3, 9;
    p17s4, p17s8, T, 4, 8;
    p18s3, p18s7, T, 3, 7;
    p19s4, p19s19, t1, 4, 19;
    p20s3, p20s19, t1, 3, 19];

%This matrix is added to in the compression scenario
Element2=Element;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Various cases specified
% Case 2 rigid body connections
CS2=[2, 0; %1
     1, 0; %2
     2, 15; %3
     1, 14; %4
     2, 17; %5
     2, 6; %6
     1, 5; %7
     1, 16; %8
     13, 6; %9
     6, 11; %10
     5, 11; %11
     12, 5; %12
      4, 6; %13
      3, 5; %14
      4, 10; %15
      3, 9; %16
      4, 8; %17
      3, 7; %18
      4, 18; %19
      3, 18]; %20
  
% Case 3 rigid body connections
CS3=[2, 0; %1
     1, 0; %2
     2, 14; %3
     1, 14; %4
     2, 16; %5
     2, 6; %6
     1, 5; %7
     1, 15; %8
     13, 6; %9
     6, 11; %10
     5, 10; %11
     12, 5; %12
      4, 6; %13
      3, 5; %14
      4, 9; %15
      3, 8; %16
      4, 7; %17
      3, 7; %18
      4, 17; %19
      3, 17]; %20
  
% Case 4 rigid body connections
CS4=[0, 0; %1
     0, 0; %2
     0, 0; %3
     0, 0; %4
     0, 1; %5
     0, 6; %6
     0, 5; %7
     0, 2; %8
     4, 6; %9
     6, 8; %10
     5, 7; %11
     3, 5; %12
     11, 6; %13
     11, 5; %14
     11, 10; %15
     11, 9; %16
     11, 11; %17
     11, 11; %18
     11, 11; %19
     11, 11]; %20
 
% Case 5 rigid body connections
CS5=[0, 0; %1
     0, 0; %2
     0, 0; %3
     0, 0; %4
     0, 1; %5
     0, 6; %6
     0, 5; %7
     0, 2; %8
     4, 6; %9
     6, 7; %10
     5, 7; %11
     3, 5; %12
     10, 6; %13
     10, 5; %14
     10, 9; %15
     10, 8; %16
     10, 10; %17
     10, 10; %18
     10, 10; %19
     10, 10]; %20

% Case 6 rigid body connections
CS6=[0, 0; %1
     0, 0; %2
     0, 0; %3
     0, 0; %4
     0, 1; %5
     0, 4; %6
     0, 3; %7
     0, 2; %8
     1, 4; %9
     4, 6; %10
     3, 5; %11
     2, 3; %12
     7, 4; %13
     7, 3; %14
     7, 1; %15
     7, 2; %16
     7, 7; %17
     7, 7; %18
     7, 7; %19
     7, 7]; %20

% Case 7 rigid body connections
CS7=[1, 0; %1
     2, 0; %2
     1, 4; %3
     2, 3; %4
     1, 5; %5
     1, 8; %6
     2, 7; %7
     2, 6; %8
     5, 8; %9
     8, 9; %10
     7, 9; %11
     6, 7; %12
     10, 8; %13
     11, 7; %14
     10, 5; %15
     11, 6; %16
     10, 12; %17
     11, 13; %18
     10, 14; %19
     11, 14]; %20 
 
% Case 8 rigid body connections
CS8=[1, 0; %1
     2, 0; %2
     1, 4; %3
     2, 3; %4
     1, 5; %5
     1, 8; %6
     2, 7; %7
     2, 6; %8
     5, 8; %9
     8, 10; %10
     7, 9; %11
     6, 7; %12
     11, 8; %13
     12, 7; %14
     11, 5; %15
     12, 6; %16
     11, 13; %17
     12, 14; %18
     11, 15; %19
     12, 15]; %20  
 
% Case 9 rigid body connections
CS9=[0, 0; %1
     0, 0; %2
     0, 0; %3
     0, 0; %4
     0, 1; %5
     0, 4; %6
     0, 3; %7
     0, 2; %8
     1, 4; %9
     4, 5; %10
     3, 5; %11
     2, 3; %12
     6, 4; %13
     6, 3; %14
     6, 1; %15
     6, 2; %16
     6, 6; %17
     6, 6; %18
     6, 6; %19
     6, 6]; %20 

% Case 10 rigid body connections
CS10=[2, 0; %1
      1, 0; %2
      2, 16; %3
      1, 15; %4
      2, 18; %5
      2, 6; %6
      1, 5; %7
      1, 17; %8
      14, 6; %9
      6, 12; %10
      5, 11; %11
      13, 5; %12
      4, 6; %13
      3, 5; %14
      4, 10; %15
      3, 9; %16
      4, 8; %17
      3, 7; %18
      4, 19; %19
      3, 19]; %20
 
% Case 11 rigid body connections
CS11=Element(1:20,8:9);
CS11(13,1:2)=[10, 6];
CS11(14,1:2)=[9, 5];
CS11(7,1:2)=[17, 5];
CS11(6,1:2)=[18, 6];

% Case 12 rigid body connections
CS12=CS2(1:20,1:2);
CS12(13,1:2)=[10, 6];
CS12(14,1:2)=[9, 5];
CS12(7,1:2)=[16, 5];
CS12(6,1:2)=[17, 6];

% Case 13 rigid body connections
CS13=CS3(1:20,1:2);
CS12(13,1:2)=[9, 6];
CS12(14,1:2)=[8, 5];
CS12(7,1:2)=[15, 5];
CS12(6,1:2)=[16, 6];

% Case 14 rigid body connections
CS14=CS4;
      
% Case 15 rigid body connections
CS15=CS5;

% Case 16 rigid body connections
CS16=CS2;

% Case 17 rigid body connections
CS17=CS10;

% Case 18 rigid body connections
CS18=CS2;

% Case 19 rigid body connections
CS19=[0, 0; %1
     0, 0; %2
     0, 0; %3
     0, 0; %4
     0, 1; %5
     0, 5; %6
     0, 4; %7
     0, 2; %8
     3, 5; %9
     5, 7; %10
     4, 6; %11
     3, 4; %12
     10, 5; %13
     10, 4; %14
     10, 9; %15
     10, 8; %16
     10, 10; %17
     10, 10; %18
     10, 10; %19
     10, 10]; %20 
 
% Case 20 rigid body connections
CS20=[0, 0; %1
     0, 0; %2
     0, 0; %3
     0, 0; %4
     0, 1; %5
     0, 5; %6
     0, 4; %7
     0, 2; %8
     3, 5; %9
     5, 6; %10
     4, 6; %11
     3, 4; %12
     9, 5; %13
     9, 4; %14
     9, 8; %15
     9, 7; %16
     9, 9; %17
     9, 9; %18
     9, 9; %19
     9, 9]; %20
 
% Case 21 rigid body connections
CS21=[0, 0; %1
     0, 0; %2
     0, 0; %3
     0, 0; %4
     0, 1; %5
     0, 3; %6
     0, 2; %7
     0, 1; %8
     1, 3; %9
     3, 5; %10
     2, 4; %11
     1, 2; %12
     6, 3; %13
     6, 2; %14
     6, 1; %15
     6, 1; %16
     6, 6; %17
     6, 6; %18
     6, 6; %19
     6, 6]; %20
 
% Case 22 rigid body connections
CS22=[0, 0; %1
     0, 0; %2
     0, 0; %3
     0, 0; %4
     0, 1; %5
     0, 3; %6
     0, 2; %7
     0, 1; %8
     1, 3; %9
     3, 4; %10
     2, 4; %11
     1, 2; %12
     5, 3; %13
     5, 2; %14
     5, 1; %15
     5, 1; %16
     5, 5; %17
     5, 5; %18
     5, 5; %19
     5, 5]; %20

% Case 23 rigid body connections
CS23=CS8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type=1; % type is one when the vertical hard stops aren't hitting the bow tie flexures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tension
g1=Element(11,4)-Element(10,4);
g2=Element(18,4)-Element(17,4);
o1=(Element(13,1:3)-Element(13,4:6))/sqrt(dot((Element(13,1:3)-Element(13,4:6)),(Element(13,1:3)-Element(13,4:6))));
o2=Element(13,4:6)-[(Element(13,7)/cos(theta2))/2, 0, 0];
o3=(Element(15,4:6)-Element(15,1:3));
no1=[-o3(1,2), o3(1,1), 0];
no2=no1/sqrt(dot(no1,no1));
o4=Element(15,4:6)+((Element(15,7)/2)*no2);
t=(o4(1,2)-o2(1,2))/o1(1,2);
g7=-(o4(1,1)-(o2(1,1)+(t*o1(1,1))));

true=1;
% Calculate the first quarter of the cylce
for j=1:n
    if(g7<=0 && true==1)
        % Correct the element thickness and locations if the vertical hardstops hit the bowtie flexures
        pp1=Element(15,((j-1)*9)+1:((j-1)*9)+3)-Element(13,((j-1)*9)+1:((j-1)*9)+3);
        pp2=pp1/sqrt(dot(pp1,pp1));
        pp3=Element(15,((j-1)*9)+1:((j-1)*9)+3)+((Element(15,((j-1)*9)+7)/2)*pp2); % left top corner of 15
        Newthick=Element(15,((j-1)*9)+7)+Element(13,((j-1)*9)+7);
        pp4=Element(15,((j-1)*9)+4:((j-1)*9)+6)-Element(15,((j-1)*9)+1:((j-1)*9)+3); % points the lentgh of 15 from top to bottom
        pp5=pp3-((Newthick/2)*pp2);
        pp6=pp5+pp4;
        pp7=Element(13,((j-1)*9)+1:((j-1)*9)+3)+pp4;
        Element(15,((j-1)*9)+1:((j-1)*9)+3)=pp5;
        Element(15,((j-1)*9)+4:((j-1)*9)+6)=pp6;
        Element(15,((j-1)*9)+7)=Newthick;
        Element(16,((j-1)*9)+1:((j-1)*9)+3)=[-pp5(1,1), pp5(1,2), 0];
        Element(16,((j-1)*9)+4:((j-1)*9)+6)=[-pp6(1,1), pp6(1,2), 0];
        Element(16,((j-1)*9)+7)=Newthick;
        Element(13,((j-1)*9)+1:((j-1)*9)+3)=pp7;
        Element(14,((j-1)*9)+1:((j-1)*9)+3)=[-pp7(1,1), pp7(1,2), 0];        
        pp8=Element(5,((j-1)*9)+1:((j-1)*9)+3)-Element(6,((j-1)*9)+1:((j-1)*9)+3);
        pp9=pp8/sqrt(dot(pp8,pp8));
        pp10=Element(5,((j-1)*9)+1:((j-1)*9)+3)+((Element(5,((j-1)*9)+7)/2)*pp9); % left bottom corner of 5
        Newthick2=Element(5,((j-1)*9)+7)+Element(6,((j-1)*9)+7);
        pp11=Element(5,((j-1)*9)+4:((j-1)*9)+6)-Element(5,((j-1)*9)+1:((j-1)*9)+3); % points the lentgh of 5 from bottom to top
        pp12=pp10-((Newthick2/2)*pp9);
        pp13=pp12+pp11;
        pp14=Element(6,((j-1)*9)+1:((j-1)*9)+3)+pp11;
        Element(5,((j-1)*9)+1:((j-1)*9)+3)=pp12;
        Element(5,((j-1)*9)+4:((j-1)*9)+6)=pp13;
        Element(5,((j-1)*9)+7)=Newthick2;
        Element(8,((j-1)*9)+1:((j-1)*9)+3)=[-pp12(1,1), pp12(1,2), 0];
        Element(8,((j-1)*9)+4:((j-1)*9)+6)=[-pp13(1,1), pp13(1,2), 0];
        Element(8,((j-1)*9)+7)=Newthick2;
        Element(6,((j-1)*9)+1:((j-1)*9)+3)=pp14;
        Element(7,((j-1)*9)+1:((j-1)*9)+3)=[-pp14(1,1), pp14(1,2), 0];
        true=2;
    end
    if(g1>0 && g2>0 && g7>0) % case 1
        type(1,j+1)=1;
        [Element,strainPoisson,g1,g2,g7]=TCase1(j,Element,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G);
    elseif(g1<=0 && g2>0 && g7>0) % case 2
        type(1,j+1)=1;
        [Element,strainPoisson,g1,g2,g7]=TCase2(j,Element,CS2,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G);
    elseif(g2<=0 && g7>0) % case 3
        type(1,j+1)=1;
        [Element,strainPoisson,g1,g2,g7]=TCase3(j,Element,CS3,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G);
    elseif(g1>0 && g2>0 && g7<=0) % case 11
        type(1,j+1)=2;
        [Element,strainPoisson,g1,g2,g7]=TCase11(j,Element,CS11,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G);
    elseif(g1<=0 && g2>0 && g7<=0) % case 12
        type(1,j+1)=2;
        [Element,strainPoisson,g1,g2,g7]=TCase12(j,Element,CS12,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G);
    elseif(g2<=0 && g7<=0) % case 13
        type(1,j+1)=2;
        [Element,strainPoisson,g1,g2,g7]=TCase13(j,Element,CS13,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G);   
    end
end
%Calculate the second quarter of the cycle
count=1;
for k=n-1:-1:0   
    strainPoisson(1,n+1+count)=strainPoisson(1,k+1);
    count=count+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compression
strainload2=-strainload;
strainPoisson2=0;

g1=Element2(11,4)-Element2(10,4);
g3=Element2(15,5)-(Element2(9,2)+(Element2(9,7)/2));
vv1=Element2(16,1:3)-Element2(14,1:3);
vv2=Element2(20,1:3)+[(Element2(20,7)/2)*sin(theta1), (Element2(20,7)/2)*cos(theta1), 0]+(((Element2(20,7)*1.5)-(Element2(20,7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element2(12,1)-vv2(1,1);
g5=Element2(19,5)+((Element2(19,7)/2)/cos(theta1))-vv2(1,2);
vspec=Element2(16,1:3)-Element2(16,4:6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element2(16,4:6)+((Element2(16,7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element2(12,1)-vv3(1,1);

% Calculate the third quarter of the cylce
for j=1:n
    if(g1>0 && g3>0 && g4>0 && g5>0 && g6>0) % case 1
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase1(j,Element2,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1<=0 && g3>0 && g4>0 && g5>0 && g6>0) % case 2
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase2(j,Element2,CS2,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1>0 && g3>0 && g4>0 && g5<=0 && g6>0) % case 4
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase4(j,Element2,CS4,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1<=0 && g3>0 && g4>0 && g5<=0 && g6>0) % case 5
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase5(j,Element2,CS5,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1>0 && g3<=0 && g4>0 && g5<=0) % case 6
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase6(j,Element2,CS6,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1<=0 && g3<=0 && g4>0 && g5>0) % case 7
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase7(j,Element2,CS7,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1>0 && g3<=0 && g4>0 && g5>0) % case 8
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase8(j,Element2,CS8,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1<=0 && g3<=0 && g4>0 && g5<=0) % case 9
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase9(j,Element2,CS9,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1>0 && g3>0 && g4>0 && g5>0 && g6<=0) % case 10
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase10(j,Element2,CS10,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH);
    elseif(g1>0 && g3>0 && g4>0 && g5<=0 && g6<=0) % case 14
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase14(j,Element2,CS14,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH);
    elseif(g1<=0 && g3>0 && g4>0 && g5<=0 && g6<=0) % case 15
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase15(j,Element2,CS15,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH);
    elseif(g1<=0 && g3>0 && g4>0 && g5>0 && g6<=0) % case 16
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase16(j,Element2,CS16,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH);
    elseif(g1>0 && g3>0 && g4<=0 && g5>0) % case 17
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase17(j,Element2,CS17,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH);
    elseif(g1<=0 && g3>0 && g4<=0 && g5>0) % case 18
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase18(j,Element2,CS18,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH);
    elseif(g1>0 && g3>0 && g4<=0 && g5<=0) % case 19
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase19(j,Element2,CS19,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1<=0 && g3>0 && g4<=0 && g5<=0) % case 20
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase20(j,Element2,CS20,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1>0 && g3<=0 && g4<=0 && g5<=0) % case 21
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase21(j,Element2,CS21,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1<=0 && g3<=0 && g4<=0 && g5<=0) % case 22
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase22(j,Element2,CS22,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n);
    elseif(g1>0 && g3<=0 && g4<=0 && g5>0) % case 23
        [Element2,strainPoisson2,g1,g3,g4,g5,g6]=CCase23(j,Element2,CS23,strainload2,strainPoisson2,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH);
    end
end
% Calculate the fourth quarter of the cycle
count=1;
for k=n-1:-1:0
    strainPoisson2(1,n+1+count)=strainPoisson2(1,k+1);
    count=count+1;
end
[r2,c2]=size(strainPoisson2);
totalstrainPoisson=[strainPoisson, strainPoisson2(1,2:c2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Draw the cell


status(2);%%Update Status
% draw all the rest of the cells in the animated gif
Output=DrawAllImages(Element,Element2,type,theta1,theta2,dim,n,S5,S6,L5,L6); %nnew
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate the Plot
figure('visible','off','color','w');
plot(time/60,100*strainload);
plot(time/60,100*totalstrainPoisson,'r');
save data time totalstrainPoisson strainload -append;
axis([0, hr/60, -1.1*la*100, 1.1*la*100]);
xlabel('Time (min)');
ylabel('Strain (%)');
Output{end+1}=getframe(gcf);
hold off

function Output=DrawAllImages(Element,Element2,type,theta1,theta2,dim,n,S5,S6,L5,L6) %nnew
% This function generates all the images that should constitute the
% animated Gif of the cell
%figure('visible','off','color','w','Unit','normalized','Position',[0 0 .5 1]);
[rrr,ccc]=size(type);
Output=cell(1,4*n);
celln=1;%Counter for output
if(type(1,ccc)==2)
    type(1,ccc+1)=2;
else
    type(1,ccc+1)=1;
end

for j=1:n
    
    %figure('visible','off','color','w','Unit','normalized','Position',[0 0 .5 1]);
    celln=celln+1;
    Output{celln}=drawcell(Element(1:20,(j*9)+1:(j*9)+9),type(1,j+2),theta1,theta2,dim,S5,S6,L5,L6); %nnew
    %Output{celln} = getframe(gca);
end

count=0;
for k=n-1:-1:0
    count=count+1;
   % figure('visible','off','color','w','Unit','normalized','Position',[0 0 .5 1]);
    celln=celln+1;
    Output{celln}=drawcell(Element(1:20,(k*9)+1:(k*9)+9),type(1,k+2),theta1,theta2,dim,S5,S6,L5,L6); %nnew
%Output{celln} = getframe(gca);
end

for j=1:n
    %figure('visible','off','color','w','Unit','normalized','Position',[0 0 .5 1]);
    celln=celln+1;
    Output{celln}=drawcell(Element2(1:20,(j*9)+1:(j*9)+9),1,theta1,theta2,dim,S5,S6,L5,L6); %nnew
%Output{celln} = getframe(gca);
end

count1=0;
for k=n-1:-1:0
    count1=count1+1;
    %figure('visible','off','color','w','Unit','normalized','Position',[0 0 .5 1]);
    celln=celln+1;
    Output{celln}=drawcell(Element2(1:20,(k*9)+1:(k*9)+9),1,theta1,theta2,dim,S5,S6,L5,L6); %nnew
%Output{celln} = getframe(gca);
end
function Output=drawcircle(x,y,r,color) %This entire function is nnew
% draws circle where x and y are center coordinates and r is the radius
% if color is 0, circle is white, if color is 1, circle is grey
format long
hold on
inc=50;
if(color==0)
    th = 0:pi/inc:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    zunit = zeros(size(yunit));
    patch(xunit, yunit, zunit, [1 1 1], 'LineStyle', 'none');
else
    th = 0:pi/inc:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    zunit = zeros(size(yunit));
    patch(xunit, yunit, zunit, [.6 .6 .6], 'LineStyle', 'none');
end
function Output=drawcell(Element,typenum,theta1,theta2,dim,S5,S6,L5,L6) %nnew
% This funtion draws the cell described by the Element matrix. typenum tells
% what type number it is. theta1 and theta2 are the paremters defined in the figure
% and dim is the dimension that sets the size of the display window

format long
fig=figure('visible','off','color','w','Position',[0 0 1600/2.3 800]);
hold on

[r,c]=size(Element);

if(typenum==1) % ordinary cases when vertical stops aren't hitting the bowtie flexures
    % Draw the elements
    for i=1:r
        drawrectangle(Element(i,1:3),Element(i,4:6)-Element(i,1:3),Element(i,7))
    end
    %%%%%%%%%%%%%%%%%% nnew
    % Draw the circular feature on element 17, 18, 3, and 4
    drawcircle(Element(17,4),Element(17,5),Element(17,7)/2,0);
    drawcircle(Element(18,4),Element(18,5),Element(18,7)/2,1);
    drawrectangle(Element(18,1:3),Element(18,4:6)-Element(18,1:3),Element(18,7))
    drawcircle(Element(3,4),Element(3,5),Element(3,7)/2,0);
    drawcircle(Element(4,4),Element(4,5),Element(4,7)/2,1);
    drawrectangle(Element(4,1:3),Element(4,4:6)-Element(4,1:3),Element(4,7))
    % Draw Tabs
    Pvec1=Element(12,1:3)+[S5/2, L5/2, 0];
    Pvec2=Element(12,1:3)+[S5/2, -L5/2, 0];
    drawrectangle(Pvec1, Pvec2-Pvec1, S5);
    Pvec3=Element(9,1:3)+[-S5/2, L5/2, 0];
    Pvec4=Element(9,1:3)+[-S5/2, -L5/2, 0];
    drawrectangle(Pvec3, Pvec4-Pvec3, S5);
    Pvec5=Element(19,4:6)+((Element(20,4:6)-Element(19,4:6))/2)+[0, Element(19,7)/(2*cos(theta1)), 0];
    Pvec6=Pvec5+[L6/2, S6/2, 0];
    Pvec7=Pvec5+[-L6/2, S6/2, 0];
    drawrectangle(Pvec6, Pvec7-Pvec6, S6);
    Pvec8=Element(1,4:6)+((Element(2,4:6)-Element(1,4:6))/2)+[0, -Element(1,7)/(2*cos(theta1)), 0];
    Pvec9=Pvec8+[L6/2, -S6/2, 0];
    Pvec10=Pvec8+[-L6/2, -S6/2, 0];
    drawrectangle(Pvec9, Pvec10-Pvec9, S6);
  %%%%%%%%%%%%%%%%%%%
    % Draw the rigid bodies
    % Draw body 0
    drawrectangle(Element(1,4:6), Element(2,4:6)-Element(1,4:6), Element(1,7)/cos(theta1));
    v1=(Element(1,4:6)-Element(1,1:3))/sqrt(dot((Element(1,4:6)-Element(1,1:3)),(Element(1,4:6)-Element(1,1:3))));
    drawrectangle(Element(1,4:6), v1*(Element(1,7)/2)*tan(theta1), Element(1,7));
    v2=[-v1(1,1), v1(1,2), 0];
    drawrectangle(Element(2,4:6), v2*(Element(2,7)/2)*tan(theta1), Element(2,7));
    % Draw body 19
    drawrectangle(Element(19,4:6), Element(20,4:6)-Element(19,4:6), Element(19,7)/cos(theta1));
    v1=[v1(1,1), -v1(1,2), 0];
    drawrectangle(Element(19,4:6), v1*(Element(19,7)/2)*tan(theta1), Element(19,7));
    v2=[-v1(1,1), v1(1,2), 0];
    drawrectangle(Element(20,4:6), v2*(Element(20,7)/2)*tan(theta1), Element(20,7));   
    % Draw body 1   
    v1=(Element(7,1:3)-Element(7,4:6))/sqrt(dot((Element(7,1:3)-Element(7,4:6)),(Element(7,1:3)-Element(7,4:6))));
    drawrectangle(Element(7,1:3), v1*(Element(7,7)/2)*tan(theta2), Element(7,7));
    v2=(Element(8,1:3)-Element(7,1:3))/sqrt(dot((Element(8,1:3)-Element(7,1:3)),(Element(8,1:3)-Element(7,1:3))));
    drawrectangle(Element(4,1:3), v2*((Element(7,7)/cos(theta2))+Element(8,7)), Element(4,7));
    drawrectangle(Element(2,1:3)+[-(Element(2,7)/2)*sin(theta1), (Element(2,7)/2)*cos(theta1), 0], v2*1.5*Element(2,7), 2*Element(2,7)*cos(theta1));  
    % Draw body 2
    v1=(Element(6,1:3)-Element(6,4:6))/sqrt(dot((Element(6,1:3)-Element(6,4:6)),(Element(6,1:3)-Element(6,4:6))));
    drawrectangle(Element(6,1:3), v1*(Element(6,7)/2)*tan(theta2), Element(6,7));
    v2=(Element(5,1:3)-Element(6,1:3))/sqrt(dot((Element(5,1:3)-Element(6,1:3)),(Element(5,1:3)-Element(6,1:3))));
    drawrectangle(Element(3,1:3), v2*((Element(6,7)/cos(theta2))+Element(5,7)), Element(3,7));
    drawrectangle(Element(1,1:3)+[(Element(1,7)/2)*sin(theta1), (Element(1,7)/2)*cos(theta1), 0], v2*1.5*Element(1,7), 2*Element(1,7)*cos(theta1));  
    % Draw body 3   
    v1=(Element(14,1:3)-Element(14,4:6))/sqrt(dot((Element(14,1:3)-Element(14,4:6)),(Element(14,1:3)-Element(14,4:6))));
    drawrectangle(Element(14,1:3), v1*(Element(14,7)/2)*tan(theta2), Element(14,7));
    v2=(Element(16,1:3)-Element(14,1:3))/sqrt(dot((Element(16,1:3)-Element(14,1:3)),(Element(16,1:3)-Element(14,1:3))));
    drawrectangle(Element(18,1:3), v2*((Element(14,7)/cos(theta2))+Element(16,7)), Element(18,7));
    drawrectangle(Element(20,1:3)+[-(Element(20,7)/2)*sin(theta1), -(Element(20,7)/2)*cos(theta1), 0], v2*1.5*Element(20,7), 2*Element(2,7)*cos(theta1));  
    % Draw body 4
    v1=(Element(13,1:3)-Element(13,4:6))/sqrt(dot((Element(13,1:3)-Element(13,4:6)),(Element(13,1:3)-Element(13,4:6))));
    drawrectangle(Element(13,1:3), v1*(Element(13,7)/2)*tan(theta2), Element(13,7));
    v2=(Element(15,1:3)-Element(13,1:3))/sqrt(dot((Element(15,1:3)-Element(13,1:3)),(Element(15,1:3)-Element(13,1:3))));
    drawrectangle(Element(17,1:3), v2*((Element(13,7)/cos(theta2))+Element(15,7)), Element(17,7));
    drawrectangle(Element(19,1:3)+[(Element(19,7)/2)*sin(theta1), -(Element(19,7)/2)*cos(theta1), 0], v2*1.5*Element(19,7), 2*Element(19,7)*cos(theta1));  
    % Draw body 5
    drawrectangle(Element(11,1:3), Element(12,4:6)-Element(11,1:3), Element(11,7));
    v1=(Element(14,4:6)-Element(14,1:3))/sqrt(dot((Element(14,4:6)-Element(14,1:3)),(Element(14,4:6)-Element(14,1:3))));
    drawrectangle(Element(14,4:6), v1*(Element(14,7)/2)*tan(theta2), Element(14,7));
    v2=(Element(7,4:6)-Element(7,1:3))/sqrt(dot((Element(7,4:6)-Element(7,1:3)),(Element(7,4:6)-Element(7,1:3))));
    drawrectangle(Element(7,4:6), v2*(Element(7,7)/2)*tan(theta2), Element(7,7));
    % Draw body 6
    drawrectangle(Element(10,1:3), Element(9,4:6)-Element(10,1:3), Element(10,7));
    v1=(Element(13,4:6)-Element(13,1:3))/sqrt(dot((Element(13,4:6)-Element(13,1:3)),(Element(13,4:6)-Element(13,1:3))));
    drawrectangle(Element(13,4:6), v1*(Element(13,7)/2)*tan(theta2), Element(13,7));
    v2=(Element(6,4:6)-Element(6,1:3))/sqrt(dot((Element(6,4:6)-Element(6,1:3)),(Element(6,4:6)-Element(6,1:3))));
    drawrectangle(Element(6,4:6), v2*(Element(6,7)/2)*tan(theta2), Element(6,7));
elseif(typenum==2) % extreme cases when vertical stops do hit the bowtie flexures
    % Draw the elements
    for i=1:r
        drawrectangle(Element(i,1:3),Element(i,4:6)-Element(i,1:3),Element(i,7))
    end
  %%%%%%%%%%%%%%%%%% nnew
    % Draw the circular feature on element 17, 18, 3, and 4
    drawcircle(Element(17,4),Element(17,5),Element(17,7)/2,0);
    drawcircle(Element(18,4),Element(18,5),Element(18,7)/2,1);
    drawcircle(Element(3,4),Element(3,5),Element(3,7)/2,0);
    drawcircle(Element(4,4),Element(4,5),Element(4,7)/2,1);
    % Draw Tabs
    Pvec1=Element(12,1:3)+[S5/2, L5/2, 0];
    Pvec2=Element(12,1:3)+[S5/2, -L5/2, 0];
    drawrectangle(Pvec1, Pvec2-Pvec1, S5);
    Pvec3=Element(9,1:3)+[-S5/2, L5/2, 0];
    Pvec4=Element(9,1:3)+[-S5/2, -L5/2, 0];
    drawrectangle(Pvec3, Pvec4-Pvec3, S5);
    Pvec5=Element(19,4:6)+((Element(20,4:6)-Element(19,4:6))/2)+[0, Element(19,7)/(2*cos(theta1)), 0];
    Pvec6=Pvec5+[L6/2, S6/2, 0];
    Pvec7=Pvec5+[-L6/2, S6/2, 0];
    drawrectangle(Pvec6, Pvec7-Pvec6, S6);
    Pvec8=Element(1,4:6)+((Element(2,4:6)-Element(1,4:6))/2)+[0, -Element(1,7)/(2*cos(theta1)), 0];
    Pvec9=Pvec8+[L6/2, -S6/2, 0];
    Pvec10=Pvec8+[-L6/2, -S6/2, 0];
    drawrectangle(Pvec9, Pvec10-Pvec9, S6);
  %%%%%%%%%%%%%%%%%%%
    % Draw the rigid bodies
    % Draw body 0
    drawrectangle(Element(1,4:6), Element(2,4:6)-Element(1,4:6), Element(1,7)/cos(theta1));
    v1=(Element(1,4:6)-Element(1,1:3))/sqrt(dot((Element(1,4:6)-Element(1,1:3)),(Element(1,4:6)-Element(1,1:3))));
    drawrectangle(Element(1,4:6), v1*(Element(1,7)/2)*tan(theta1), Element(1,7));
    v2=[-v1(1,1), v1(1,2), 0];
    drawrectangle(Element(2,4:6), v2*(Element(2,7)/2)*tan(theta1), Element(2,7));
    % Draw body 19
    drawrectangle(Element(19,4:6), Element(20,4:6)-Element(19,4:6), Element(19,7)/cos(theta1));
    v1=[v1(1,1), -v1(1,2), 0];
    drawrectangle(Element(19,4:6), v1*(Element(19,7)/2)*tan(theta1), Element(19,7));
    v2=[-v1(1,1), v1(1,2), 0];
    drawrectangle(Element(20,4:6), v2*(Element(20,7)/2)*tan(theta1), Element(20,7));   
    % Draw body 4
    aa1=[Element(15,1)-Element(17,1); Element(15,2)-Element(17,2)];
    mat1=[Element(17,7)/2, -Element(15,7)/2; -Element(15,7)/2, -Element(17,7)/2];  
    bb1=mat1\aa1;
    angle1=asin(bb1(1,1));
    v1=-[cos(angle1), sin(angle1), 0];
    drawrectangle(Element(17,1:3), Element(15,7)*v1, Element(17,7));
    drawrectangle(Element(19,1:3)+[(Element(19,7)/2)*sin(theta1), -(Element(19,7)/2)*cos(theta1), 0], v1*1.5*Element(19,7), 2*Element(19,7)*cos(theta1));     
    % Draw body 1   
    v2=[-v1(1,1), -v1(1,2), 0]; 
    drawrectangle(Element(4,1:3), Element(8,7)*v2, Element(4,7));
    drawrectangle(Element(2,1:3)+[-(Element(2,7)/2)*sin(theta1), (Element(2,7)/2)*cos(theta1), 0], v2*1.5*Element(2,7), 2*Element(2,7)*cos(theta1));  
    % Draw body 2
    v2=[v1(1,1),-v1(1,2),0];
    drawrectangle(Element(3,1:3), Element(5,7)*v2, Element(3,7));
    drawrectangle(Element(1,1:3)+[(Element(1,7)/2)*sin(theta1), (Element(1,7)/2)*cos(theta1), 0], v2*1.5*Element(1,7), 2*Element(1,7)*cos(theta1));  
    % Draw body 3   
    v2=[-v1(1,1),v1(1,2),0];
    drawrectangle(Element(18,1:3), Element(16,7)*v2, Element(18,7));
    drawrectangle(Element(20,1:3)+[-(Element(20,7)/2)*sin(theta1), -(Element(20,7)/2)*cos(theta1), 0], v2*1.5*Element(20,7), 2*Element(20,7)*cos(theta1));
    % Draw body 5
    drawrectangle(Element(11,1:3), Element(12,4:6)-Element(11,1:3), Element(11,7));
    v1=(Element(14,4:6)-Element(14,1:3))/sqrt(dot((Element(14,4:6)-Element(14,1:3)),(Element(14,4:6)-Element(14,1:3))));
    drawrectangle(Element(14,4:6), v1*(Element(14,7)/2)*tan(theta2), Element(14,7));
    v2=(Element(7,4:6)-Element(7,1:3))/sqrt(dot((Element(7,4:6)-Element(7,1:3)),(Element(7,4:6)-Element(7,1:3))));
    drawrectangle(Element(7,4:6), v2*(Element(7,7)/2)*tan(theta2), Element(7,7));
    % Draw body 6
    drawrectangle(Element(10,1:3), Element(9,4:6)-Element(10,1:3), Element(10,7));
    v1=(Element(13,4:6)-Element(13,1:3))/sqrt(dot((Element(13,4:6)-Element(13,1:3)),(Element(13,4:6)-Element(13,1:3))));
    drawrectangle(Element(13,4:6), v1*(Element(13,7)/2)*tan(theta2), Element(13,7));
    v2=(Element(6,4:6)-Element(6,1:3))/sqrt(dot((Element(6,4:6)-Element(6,1:3)),(Element(6,4:6)-Element(6,1:3))));
    drawrectangle(Element(6,4:6), v2*(Element(6,7)/2)*tan(theta2), Element(6,7));  
end
hold off
axis on
axis equal
axis([-dim, dim, -dim, dim*1.3]);
xlabel('x-axis (m)');
ylabel('y-axis (m)');
set(gca,'FontSize',25);
Output = getframe(fig);
function drawrectangle(Loc,Lsvec,t)
% This funtion draws a rectangular strap of thickness t that begins at Loc 
% and points in the direction and length of Lsvec. It also adds an extra 
% length to both sides of length extra.

format long

hold on
extra=0;

Ls=sqrt(dot(Lsvec,Lsvec));
unitLs=Lsvec/Ls;
perpLs=[unitLs(1,2), -unitLs(1,1)];
strapx=[(Loc(1,1)+((t/2)*perpLs(1,1))-(extra*unitLs(1,1))), (Loc(1,1)-((t/2)*perpLs(1,1))-(extra*unitLs(1,1))), (Loc(1,1)-((t/2)*perpLs(1,1))+((Ls+extra)*unitLs(1,1))), (Loc(1,1)+((t/2)*perpLs(1,1))+((Ls+extra)*unitLs(1,1)))]; 
strapy=[(Loc(1,2)+((t/2)*perpLs(1,2))-(extra*unitLs(1,2))), (Loc(1,2)-((t/2)*perpLs(1,2))-(extra*unitLs(1,2))), (Loc(1,2)-((t/2)*perpLs(1,2))+((Ls+extra)*unitLs(1,2))), (Loc(1,2)+((t/2)*perpLs(1,2))+((Ls+extra)*unitLs(1,2)))]; 
strapz = zeros(size(strapy));


patch(strapx, strapy, strapz, [0.6 0.6 0.6], 'LineStyle', 'none');

hold off
function [KK]=EulerStiffnessMatrixConnect(Constraint1)
%Function reads in a matrix (Constraint) that contains all the info needed
%to describe a flexure system, and returns its stiffness matrix

%Before you construct the matrix Constraint, be sure to number all of your
%system's rigid stages and number all of your flexure elements.  Grounded 
%stages are numbered zero.

%Each row of the Constraint matrix corresponds to a flexure element.

%The following describes what each component of an individual row entails.
%the first 3 components are the chosen location vector, L, that points to
%one end of the element (the end doesn't matter but be consistent)
%the next 3 components are the direction of the n3 unit vector that points
%into the stage that the element attaches to along the element's axis.
%the next 3 components are the direction of the n2 unit vector that points
%perpendicular to the width, W, side of the rectangular element.  W is the
%length of the flat face of a blade flexure.
%the 10th component is the length of the element, L.
%the 11th component is the width of the element, W.
%the 12th component is the thickness of the element, t.
%the 13th component is the Young's Modulus of the material, E.
%the 14th component is the Shear Modulus of the material, G.
%the last two components correspond to what stages the element joings.
%start with the stage number that the chosen L vector points to and then
%end with the other stage that the element joins.

%Once you have filled a row with every element in the system make a final
%row that contains the following components:
%the 1st component is the number of stages in the system (excluding the
%ground, which is numbered 0).
%the 2nd component will always be zero.
%the 3rd component is how many rectangular elements are in the entire 
%system.  Then fill the rest of the row with zeros.


format long;

[row, column]=size(Constraint1);
stages=Constraint1(row,1);
wire=Constraint1(row,2);
blade=Constraint1(row,3);
ConstNum=row-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% take out constraints that are connected to the same stage
history=0;
countt=0;
for gg=1:ConstNum
    if(Constraint1(gg,column-1)==Constraint1(gg,column))
        countt=countt+1;
        history(1,countt)=gg;
    end
end
for hh=1:countt
    if(history(1,hh)<=wire)
        wire=wire-1;
    end
end
for ii=1:countt
    if(history(1,ii)>wire)
        blade=blade-1;
    end
end
FF=zeros(wire+blade+1,column);
counttt=0;
for gg=1:ConstNum
    if(Constraint1(gg,column-1)~=Constraint1(gg,column))
        counttt=counttt+1;
        FF(counttt,1:column)=Constraint1(gg,1:column);
    end
end
Constraint=FF;
Constraint(wire+blade+1,1)=stages;
Constraint(wire+blade+1,2)=wire;
Constraint(wire+blade+1,3)=blade;
[row, column]=size(Constraint);
ConstNum=row-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MainMatrix=zeros(stages*6,stages*6);

NumStageConst=zeros(stages,1);
for i=1:ConstNum
    for j=1:stages
        if(Constraint(i,column-1)==j || Constraint(i,column)==j)
            NumStageConst(j,1)=NumStageConst(j,1)+1;
        end
    end
end

for j=1:stages
    Iden=zeros(6,6*NumStageConst(j,1));
    for k=0:NumStageConst(j,1)-1
        Iden(1:6,(k*6)+1:(k*6)+6)=eye(6);
    end
    WrenchMatrix=zeros(6*NumStageConst(j,1),6*stages);
    count=0;
    for i=1:wire
        if(Constraint(i,column-1)==j)
            count=count+1;
            %Construct S for wire flexure i
            ss=zeros(6,6);
            l=Constraint(i,10);
            d=Constraint(i,11);
            E=Constraint(i,13);
            G=Constraint(i,14);
            I=(pi*(d^4))/64;
            J=(pi*(d^4))/32;
            A=(pi*(d^2))/4;
            ss(1,1)=l/(E*I);
            ss(1,5)=-(l^2)/(2*E*I);
            ss(2,2)=l/(E*I);
            ss(2,4)=(l^2)/(2*E*I);
            ss(3,3)=l/(G*J);
            ss(4,2)=(l^2)/(2*E*I);
            ss(4,4)=(l^3)/(3*E*I);
            ss(5,1)=-(l^2)/(2*E*I);
            ss(5,5)=(l^3)/(3*E*I);
            ss(6,6)=l/(E*A);
            S=inv(ss);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            r=Constraint(i,1:3);
            n3=Constraint(i,4:6);
            n3=n3/sqrt(dot(n3,n3));
            orth=null(n3);
            n2=transpose(orth(1:3,1));
            n2=n2/sqrt(dot(n2,n2));
            n1=cross(n2,n3);
            Na=zeros(6,6);
            Na(1:3,1)=transpose(n1);
            Na(1:3,2)=transpose(n2);
            Na(1:3,3)=transpose(n3);
            Na(4:6,4)=transpose(n1);
            Na(4:6,5)=transpose(n2);
            Na(4:6,6)=transpose(n3);
            Na(4:6,1)=transpose(cross(r,n1));
            Na(4:6,2)=transpose(cross(r,n2));
            Na(4:6,3)=transpose(cross(r,n3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            NR=zeros(6,6);
            NR(1:6,1:3)=Na(1:6,4:6);
            NR(1:6,4:6)=Na(1:6,1:3);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Nb=zeros(6,6);
            Nb(1:3,1)=transpose(n1);
            Nb(1:3,2)=transpose(n2);
            Nb(1:3,3)=transpose(n3);
            Nb(4:6,4)=transpose(n1);
            Nb(4:6,5)=transpose(n2);
            Nb(4:6,6)=transpose(n3);
            Nb(4:6,1)=transpose(cross((r-(l*n3)),n1));
            Nb(4:6,2)=transpose(cross((r-(l*n3)),n2));
            Nb(4:6,3)=transpose(cross((r-(l*n3)),n3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            WrenchMatrix(((count-1)*6)+1:((count-1)*6)+6,((j-1)*6)+1:((j-1)*6)+6)=NR*S/(Na);
            if(Constraint(i,column)~=0)
                P=zeros(6,6);
                P(4,2)=-l;
                P(5,1)=l;
                WrenchMatrix(((count-1)*6)+1:((count-1)*6)+6,((Constraint(i,column)-1)*6)+1:((Constraint(i,column)-1)*6)+6)=NR*S*(P-eye(6))/(Nb);
            end
        elseif(Constraint(i,column)==j)
            count=count+1;
            %Construct S for wire flexure i
            ss=zeros(6,6);
            l=Constraint(i,10);
            d=Constraint(i,11);
            E=Constraint(i,13);
            G=Constraint(i,14);
            I=(pi*(d^4))/64;
            J=(pi*(d^4))/32;
            A=(pi*(d^2))/4;
            ss(1,1)=l/(E*I);
            ss(1,5)=-(l^2)/(2*E*I);
            ss(2,2)=l/(E*I);
            ss(2,4)=(l^2)/(2*E*I);
            ss(3,3)=l/(G*J);
            ss(4,2)=(l^2)/(2*E*I);
            ss(4,4)=(l^3)/(3*E*I);
            ss(5,1)=-(l^2)/(2*E*I);
            ss(5,5)=(l^3)/(3*E*I);
            ss(6,6)=l/(E*A);
            S=inv(ss);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n3=-Constraint(i,4:6);
            n3=n3/sqrt(dot(n3,n3));
            r=Constraint(i,1:3)+l*n3;
            orth=null(n3);
            n2=transpose(orth(1:3,1));
            n2=n2/sqrt(dot(n2,n2));
            n1=cross(n2,n3);
            Na=zeros(6,6);
            Na(1:3,1)=transpose(n1);
            Na(1:3,2)=transpose(n2);
            Na(1:3,3)=transpose(n3);
            Na(4:6,4)=transpose(n1);
            Na(4:6,5)=transpose(n2);
            Na(4:6,6)=transpose(n3);
            Na(4:6,1)=transpose(cross(r,n1));
            Na(4:6,2)=transpose(cross(r,n2));
            Na(4:6,3)=transpose(cross(r,n3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            NR=zeros(6,6);
            NR(1:6,1:3)=Na(1:6,4:6);
            NR(1:6,4:6)=Na(1:6,1:3);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Nb=zeros(6,6);
            Nb(1:3,1)=transpose(n1);
            Nb(1:3,2)=transpose(n2);
            Nb(1:3,3)=transpose(n3);
            Nb(4:6,4)=transpose(n1);
            Nb(4:6,5)=transpose(n2);
            Nb(4:6,6)=transpose(n3);
            Nb(4:6,1)=transpose(cross((r-(l*n3)),n1));
            Nb(4:6,2)=transpose(cross((r-(l*n3)),n2));
            Nb(4:6,3)=transpose(cross((r-(l*n3)),n3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            WrenchMatrix(((count-1)*6)+1:((count-1)*6)+6,((j-1)*6)+1:((j-1)*6)+6)=NR*S/(Na);
            if(Constraint(i,column-1)~=0)
                P=zeros(6,6);
                P(4,2)=-l;
                P(5,1)=l;
                WrenchMatrix(((count-1)*6)+1:((count-1)*6)+6,((Constraint(i,column-1)-1)*6)+1:((Constraint(i,column-1)-1)*6)+6)=NR*S*(P-eye(6))/(Nb);
            end
        end
    end
    for i=1+wire:blade+wire
        if(Constraint(i,column-1)==j)
            count=count+1;
            %Construct S for blade flexure i
            ss=zeros(6,6);
            l=Constraint(i,10);
            w=Constraint(i,11);
            t=Constraint(i,12);
            E=Constraint(i,13);
            G=Constraint(i,14);
            I1=w*(t^3)/12;
            I2=t*(w^3)/12;
            Temp=0;
            if(w>t)
                for n=1:2:7
                    Temp=Temp+(tanh(n*pi*w/(2*t))/(n^5));
                end
                J=((t^3)*w/3)*(1-((192*t/((pi^5)*w))*Temp));
            else
                for n=1:2:7
                    Temp=Temp+(tanh(n*pi*t/(2*w))/(n^5));
                end
                J=((w^3)*t/3)*(1-((192*w/((pi^5)*t))*Temp));
            end
            A=w*t;
            ss(1,1)=l/(E*I1);
            ss(1,5)=-(l^2)/(2*E*I1);
            ss(2,2)=l/(E*I2);
            ss(2,4)=(l^2)/(2*E*I2);
            ss(3,3)=l/(G*J);
            ss(4,2)=(l^2)/(2*E*I2);
            ss(4,4)=(l^3)/(3*E*I2);
            ss(5,1)=-(l^2)/(2*E*I1);
            ss(5,5)=(l^3)/(3*E*I1);
            ss(6,6)=l/(E*A);
            S=inv(ss);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            r=Constraint(i,1:3);
            n3=Constraint(i,4:6);
            n3=n3/sqrt(dot(n3,n3));
            n2=Constraint(i,7:9);
            n2=n2/sqrt(dot(n2,n2));
            n1=cross(n2,n3);
            Na=zeros(6,6);
            Na(1:3,1)=transpose(n1);
            Na(1:3,2)=transpose(n2);
            Na(1:3,3)=transpose(n3);
            Na(4:6,4)=transpose(n1);
            Na(4:6,5)=transpose(n2);
            Na(4:6,6)=transpose(n3);
            Na(4:6,1)=transpose(cross(r,n1));
            Na(4:6,2)=transpose(cross(r,n2));
            Na(4:6,3)=transpose(cross(r,n3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            NR=zeros(6,6);
            NR(1:6,1:3)=Na(1:6,4:6);
            NR(1:6,4:6)=Na(1:6,1:3);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Nb=zeros(6,6);
            Nb(1:3,1)=transpose(n1);
            Nb(1:3,2)=transpose(n2);
            Nb(1:3,3)=transpose(n3);
            Nb(4:6,4)=transpose(n1);
            Nb(4:6,5)=transpose(n2);
            Nb(4:6,6)=transpose(n3);
            Nb(4:6,1)=transpose(cross((r-(l*n3)),n1));
            Nb(4:6,2)=transpose(cross((r-(l*n3)),n2));
            Nb(4:6,3)=transpose(cross((r-(l*n3)),n3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            WrenchMatrix(((count-1)*6)+1:((count-1)*6)+6,((j-1)*6)+1:((j-1)*6)+6)=NR*S/(Na);
            if(Constraint(i,column)~=0)
                P=zeros(6,6);
                P(4,2)=-l;
                P(5,1)=l;
                WrenchMatrix(((count-1)*6)+1:((count-1)*6)+6,((Constraint(i,column)-1)*6)+1:((Constraint(i,column)-1)*6)+6)=NR*S*(P-eye(6))/(Nb);
            end
        elseif(Constraint(i,column)==j)
            count=count+1;
            %Construct S for blade flexure i
            ss=zeros(6,6);
            l=Constraint(i,10);
            w=Constraint(i,11);
            t=Constraint(i,12);
            E=Constraint(i,13);
            G=Constraint(i,14);
            I1=w*(t^3)/12;
            I2=t*(w^3)/12;
            Temp=0;
            if(w>t)
                for n=1:2:7
                    Temp=Temp+(tanh(n*pi*w/(2*t))/(n^5));
                end
                J=((t^3)*w/3)*(1-((192*t/((pi^5)*w))*Temp));
            else
                for n=1:2:7
                    Temp=Temp+(tanh(n*pi*t/(2*w))/(n^5));
                end
                J=((w^3)*t/3)*(1-((192*w/((pi^5)*t))*Temp));
            end
            A=w*t;
            ss(1,1)=l/(E*I1);
            ss(1,5)=-(l^2)/(2*E*I1);
            ss(2,2)=l/(E*I2);
            ss(2,4)=(l^2)/(2*E*I2);
            ss(3,3)=l/(G*J);
            ss(4,2)=(l^2)/(2*E*I2);
            ss(4,4)=(l^3)/(3*E*I2);
            ss(5,1)=-(l^2)/(2*E*I1);
            ss(5,5)=(l^3)/(3*E*I1);
            ss(6,6)=l/(E*A);
            S=inv(ss);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n3=-Constraint(i,4:6);
            n3=n3/sqrt(dot(n3,n3));
            r=Constraint(i,1:3)+l*n3;
            n2=Constraint(i,7:9);
            n2=n2/sqrt(dot(n2,n2));
            n1=cross(n2,n3);
            Na=zeros(6,6);
            Na(1:3,1)=transpose(n1);
            Na(1:3,2)=transpose(n2);
            Na(1:3,3)=transpose(n3);
            Na(4:6,4)=transpose(n1);
            Na(4:6,5)=transpose(n2);
            Na(4:6,6)=transpose(n3);
            Na(4:6,1)=transpose(cross(r,n1));
            Na(4:6,2)=transpose(cross(r,n2));
            Na(4:6,3)=transpose(cross(r,n3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            NR=zeros(6,6);
            NR(1:6,1:3)=Na(1:6,4:6);
            NR(1:6,4:6)=Na(1:6,1:3);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Nb=zeros(6,6);
            Nb(1:3,1)=transpose(n1);
            Nb(1:3,2)=transpose(n2);
            Nb(1:3,3)=transpose(n3);
            Nb(4:6,4)=transpose(n1);
            Nb(4:6,5)=transpose(n2);
            Nb(4:6,6)=transpose(n3);
            Nb(4:6,1)=transpose(cross((r-(l*n3)),n1));
            Nb(4:6,2)=transpose(cross((r-(l*n3)),n2));
            Nb(4:6,3)=transpose(cross((r-(l*n3)),n3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            WrenchMatrix(((count-1)*6)+1:((count-1)*6)+6,((j-1)*6)+1:((j-1)*6)+6)=NR*S/(Na);
            if(Constraint(i,column-1)~=0)
                P=zeros(6,6);
                P(4,2)=-l;
                P(5,1)=l;
                WrenchMatrix(((count-1)*6)+1:((count-1)*6)+6,((Constraint(i,column-1)-1)*6)+1:((Constraint(i,column-1)-1)*6)+6)=NR*S*(P-eye(6))/(Nb);
            end            
        end
    end
    MainMatrix(((j-1)*6)+1:((j-1)*6)+6,1:(stages*6))=Iden*WrenchMatrix;
end

KK=MainMatrix;

function [Element,strainPoisson,g1,g2,g7]=TCase1(j,Element,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=19;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=[Element(i,((j-1)*9)+8), Element(i,((j-1)*9)+9)];
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=E*la*(dim^2);
max=Force;
min=0;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((max-min)/2)+min; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp>Dispinc)
        Force=((max-min)/2)+min;
        max=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    elseif(Disp<Dispinc)
        Force=((max-min)/2)+min;
        min=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    else
        Force=((max-min)/2)+min;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];
% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(Element(i,((j-1)*9)+8)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*Element(i,((j-1)*9)+8))-6+1:(6*Element(i,((j-1)*9)+8)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(Element(i,((j-1)*9)+9)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*Element(i,((j-1)*9)+9))-6+1:(6*Element(i,((j-1)*9)+9)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7:(j*9)+9)=Element(1:20,((j-1)*9)+7:((j-1)*9)+9);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g2=Element(18,(j*9)+4)-Element(17,(j*9)+4);
o1=(Element(13,(j*9)+1:(j*9)+3)-Element(13,(j*9)+4:(j*9)+6))/sqrt(dot((Element(13,(j*9)+1:(j*9)+3)-Element(13,(j*9)+4:(j*9)+6)),(Element(13,(j*9)+1:(j*9)+3)-Element(13,(j*9)+4:(j*9)+6))));
o2=Element(13,(j*9)+4:(j*9)+6)-[(Element(13,(j*9)+7)/cos(theta2))/2, 0, 0];
o3=(Element(15,(j*9)+4:(j*9)+6)-Element(15,(j*9)+1:(j*9)+3));
no1=[-o3(1,2), o3(1,1), 0];
no2=no1/sqrt(dot(no1,no1));
o4=Element(15,(j*9)+4:(j*9)+6)+((Element(15,(j*9)+7)/2)*no2);
t=(o4(1,2)-o2(1,2))/o1(1,2);
g7=-(o4(1,1)-(o2(1,1)+(t*o1(1,1))));

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g2,g7]=TCase2(j,Element,CS2,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=18;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS2(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=E*la*(dim^2);
max=Force;
min=0;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((max-min)/2)+min; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp>Dispinc)
        Force=((max-min)/2)+min;
        max=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    elseif(Disp<Dispinc)
        Force=((max-min)/2)+min;
        min=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    else
        Force=((max-min)/2)+min;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];
% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(Element(i,((j-1)*9)+8)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS2(i,1))-6+1:(6*CS2(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(Element(i,((j-1)*9)+9)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS2(i,2))-6+1:(6*CS2(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS2(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g2=Element(18,(j*9)+4)-Element(17,(j*9)+4);
o1=(Element(13,(j*9)+1:(j*9)+3)-Element(13,(j*9)+4:(j*9)+6))/sqrt(dot((Element(13,(j*9)+1:(j*9)+3)-Element(13,(j*9)+4:(j*9)+6)),(Element(13,(j*9)+1:(j*9)+3)-Element(13,(j*9)+4:(j*9)+6))));
o2=Element(13,(j*9)+4:(j*9)+6)-[(Element(13,(j*9)+7)/cos(theta2))/2, 0, 0];
o3=(Element(15,(j*9)+4:(j*9)+6)-Element(15,(j*9)+1:(j*9)+3));
no1=[-o3(1,2), o3(1,1), 0];
no2=no1/sqrt(dot(no1,no1));
o4=Element(15,(j*9)+4:(j*9)+6)+((Element(15,(j*9)+7)/2)*no2);
t=(o4(1,2)-o2(1,2))/o1(1,2);
g7=-(o4(1,1)-(o2(1,1)+(t*o1(1,1))));

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g2,g7]=TCase3(j,Element,CS3,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=17;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS3(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=E*la*(dim^2);
max=Force;
min=0;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((max-min)/2)+min; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp>Dispinc)
        Force=((max-min)/2)+min;
        max=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    elseif(Disp<Dispinc)
        Force=((max-min)/2)+min;
        min=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    else
        Force=((max-min)/2)+min;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];
% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(Element(i,((j-1)*9)+8)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS3(i,1))-6+1:(6*CS3(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(Element(i,((j-1)*9)+9)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS3(i,2))-6+1:(6*CS3(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS3(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g2=Element(18,(j*9)+4)-Element(17,(j*9)+4);
o1=(Element(13,(j*9)+1:(j*9)+3)-Element(13,(j*9)+4:(j*9)+6))/sqrt(dot((Element(13,(j*9)+1:(j*9)+3)-Element(13,(j*9)+4:(j*9)+6)),(Element(13,(j*9)+1:(j*9)+3)-Element(13,(j*9)+4:(j*9)+6))));
o2=Element(13,(j*9)+4:(j*9)+6)-[(Element(13,(j*9)+7)/cos(theta2))/2, 0, 0];
o3=(Element(15,(j*9)+4:(j*9)+6)-Element(15,(j*9)+1:(j*9)+3));
no1=[-o3(1,2), o3(1,1), 0];
no2=no1/sqrt(dot(no1,no1));
o4=Element(15,(j*9)+4:(j*9)+6)+((Element(15,(j*9)+7)/2)*no2);
t=(o4(1,2)-o2(1,2))/o1(1,2);
g7=-(o4(1,1)-(o2(1,1)+(t*o1(1,1))));

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g2,g7]=TCase11(j,Element,CS11,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=19;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS11(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=E*la*(dim^2);
max=Force;
min=0;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((max-min)/2)+min; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp>Dispinc)
        Force=((max-min)/2)+min;
        max=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    elseif(Disp<Dispinc)
        Force=((max-min)/2)+min;
        min=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    else
        Force=((max-min)/2)+min;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];
% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(Element(i,((j-1)*9)+8)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS11(i,1))-6+1:(6*CS11(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(Element(i,((j-1)*9)+9)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS11(i,2))-6+1:(6*CS11(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS11(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g2=Element(18,(j*9)+4)-Element(17,(j*9)+4);
g7=0;

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g2,g7]=TCase12(j,Element,CS12,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=18;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS12(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=E*la*(dim^2);
max=Force;
min=0;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((max-min)/2)+min; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp>Dispinc)
        Force=((max-min)/2)+min;
        max=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    elseif(Disp<Dispinc)
        Force=((max-min)/2)+min;
        min=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    else
        Force=((max-min)/2)+min;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];
% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(Element(i,((j-1)*9)+8)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS12(i,1))-6+1:(6*CS12(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(Element(i,((j-1)*9)+9)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS12(i,2))-6+1:(6*CS12(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS12(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g2=Element(18,(j*9)+4)-Element(17,(j*9)+4);
g7=0;

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g2,g7]=TCase13(j,Element,CS13,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=17;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS13(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=E*la*(dim^2);
max=Force;
min=0;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((max-min)/2)+min; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp>Dispinc)
        Force=((max-min)/2)+min;
        max=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    elseif(Disp<Dispinc)
        Force=((max-min)/2)+min;
        min=((max-min)/2)+min;
        error=abs(Disp-Dispinc);
    else
        Force=((max-min)/2)+min;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];
% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(Element(i,((j-1)*9)+8)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS13(i,1))-6+1:(6*CS13(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(Element(i,((j-1)*9)+9)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS13(i,2))-6+1:(6*CS13(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS13(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g2=Element(18,(j*9)+4)-Element(17,(j*9)+4);
g7=0;

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase1(j,Element,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,~)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=19;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=[Element(i,((j-1)*9)+8), Element(i,((j-1)*9)+9)];
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force1; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(Element(i,((j-1)*9)+8)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*Element(i,((j-1)*9)+8))-6+1:(6*Element(i,((j-1)*9)+8)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(Element(i,((j-1)*9)+9)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*Element(i,((j-1)*9)+9))-6+1:(6*Element(i,((j-1)*9)+9)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7:(j*9)+9)=Element(1:20,((j-1)*9)+7:((j-1)*9)+9);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase2(j,Element,CS2,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=18;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS2(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force1; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS2(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS2(i,1))-6+1:(6*CS2(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS2(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS2(i,2))-6+1:(6*CS2(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS2(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);


% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase4(j,Element,CS4,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=11;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS4(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force1; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS4(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS4(i,1))-6+1:(6*CS4(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS4(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS4(i,2))-6+1:(6*CS4(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS4(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase5(j,Element,CS5,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=10;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS5(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force1; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS5(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS5(i,1))-6+1:(6*CS5(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS5(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS5(i,2))-6+1:(6*CS5(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS5(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase6(j,Element,CS6,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=7;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS6(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS6(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS6(i,1))-6+1:(6*CS6(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS6(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS6(i,2))-6+1:(6*CS6(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS6(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase7(j,Element,CS7,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=14;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS7(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS7(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS7(i,1))-6+1:(6*CS7(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS7(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS7(i,2))-6+1:(6*CS7(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS7(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase8(j,Element,CS8,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=15;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS8(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS8(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS8(i,1))-6+1:(6*CS8(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS8(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS8(i,2))-6+1:(6*CS8(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS8(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase9(j,Element,CS9,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=6;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS9(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS9(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS9(i,1))-6+1:(6*CS9(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS9(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS9(i,2))-6+1:(6*CS9(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS9(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase10(j,Element,CS10,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=19;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS10(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

bb1=19;
bb2=13;
bb3=14;
bb4=9;
bb5=17;
bb6=10;
bb7=18;
vspec=Element(16,((j-1)*9)+1:((j-1)*9)+3)-Element(16,((j-1)*9)+4:((j-1)*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vector1=Element(16,((j-1)*9)+4:((j-1)*9)+6)+((Element(16,((j-1)*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
vspec=Element(8,((j-1)*9)+1:((j-1)*9)+3)-Element(8,((j-1)*9)+4:((j-1)*9)+6);
vspecn=[-vspec(1,2), vspec(1,1), 0];
vector2=Element(8,((j-1)*9)+4:((j-1)*9)+6)+((Element(8,((j-1)*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fvv=[((min-max)/2)+max, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    N=eye(6);
    N(4:6,1)=transpose(cross(vector1,[1 0 0]));
    N(4:6,2)=transpose(cross(vector1,[0 1 0]));
    N(4:6,3)=transpose(cross(vector1,[0 0 1]));
    NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
    Disp1=NewT(4,1);
    N=eye(6);
    N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
    N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
    N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
    NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
    Disp2=-NewT(4,1);
    Disp=Disp1+Disp2;
    Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
    if(Disp<Dispinc)
        Force2=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force2=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force2=((min-max)/2)+max;
        error=0;
    end
end
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fv=[0, ((min-max)/2)+max, 0];
    Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
    Fvv=[Force2, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    Disp=Twist(((6*bb1)-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
for hhh=1:HH
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, Force1, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[((min-max)/2)+max, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        N=eye(6);
        N(4:6,1)=transpose(cross(vector1,[1 0 0]));
        N(4:6,2)=transpose(cross(vector1,[0 1 0]));
        N(4:6,3)=transpose(cross(vector1,[0 0 1]));
        NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
        Disp1=NewT(4,1);
        N=eye(6);
        N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
        Disp2=-NewT(4,1);
        Disp=Disp1+Disp2;
        Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
        if(Disp<Dispinc)
            Force2=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force2=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force2=((min-max)/2)+max;
            error=0;
        end
    end
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, ((min-max)/2)+max, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[Force2, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        Disp=Twist(((6*bb1)-6)+5,1);
        Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
        if(Disp<Dispinc)
            Force1=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force1=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force1=((min-max)/2)+max;
            error=0;
        end
    end
end
Wrench=zeros(6*Constraint(21,1),1);
Fv=[0, Force1, 0];
Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
Fvv=[Force2, 0, 0];
Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS10(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS10(i,1))-6+1:(6*CS10(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS10(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS10(i,2))-6+1:(6*CS10(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS10(1:20,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase14(j,Element,CS14,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=11;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS14(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

bb1=11;
bb2=3;
bb3=4;
bb4=9;
bb5=2;
bb6=10;
bb7=1;
vspec=Element(16,((j-1)*9)+1:((j-1)*9)+3)-Element(16,((j-1)*9)+4:((j-1)*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vector1=Element(16,((j-1)*9)+4:((j-1)*9)+6)+((Element(16,((j-1)*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
vspec=Element(8,((j-1)*9)+1:((j-1)*9)+3)-Element(8,((j-1)*9)+4:((j-1)*9)+6);
vspecn=[-vspec(1,2), vspec(1,1), 0];
vector2=Element(8,((j-1)*9)+4:((j-1)*9)+6)+((Element(8,((j-1)*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fvv=[((min-max)/2)+max, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    N=eye(6);
    N(4:6,1)=transpose(cross(vector1,[1 0 0]));
    N(4:6,2)=transpose(cross(vector1,[0 1 0]));
    N(4:6,3)=transpose(cross(vector1,[0 0 1]));
    NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
    Disp1=NewT(4,1);
    N=eye(6);
    N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
    N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
    N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
    NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
    Disp2=-NewT(4,1);
    Disp=Disp1+Disp2;
    Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
    if(Disp<Dispinc)
        Force2=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force2=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force2=((min-max)/2)+max;
        error=0;
    end
end
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fv=[0, ((min-max)/2)+max, 0];
    Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
    Fvv=[Force2, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    Disp=Twist(((6*bb1)-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
for hhh=1:HH
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, Force1, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[((min-max)/2)+max, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        N=eye(6);
        N(4:6,1)=transpose(cross(vector1,[1 0 0]));
        N(4:6,2)=transpose(cross(vector1,[0 1 0]));
        N(4:6,3)=transpose(cross(vector1,[0 0 1]));
        NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
        Disp1=NewT(4,1);
        N=eye(6);
        N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
        Disp2=-NewT(4,1);
        Disp=Disp1+Disp2;
        Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
        if(Disp<Dispinc)
            Force2=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force2=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force2=((min-max)/2)+max;
            error=0;
        end
    end
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, ((min-max)/2)+max, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[Force2, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        Disp=Twist(((6*bb1)-6)+5,1);
        Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
        if(Disp<Dispinc)
            Force1=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force1=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force1=((min-max)/2)+max;
            error=0;
        end
    end
end
Wrench=zeros(6*Constraint(21,1),1);
Fv=[0, Force1, 0];
Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
Fvv=[Force2, 0, 0];
Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS14(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS14(i,1))-6+1:(6*CS14(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS14(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS14(i,2))-6+1:(6*CS14(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS14(1:20,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase15(j,Element,CS15,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=10;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS15(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

bb1=10;
bb2=3;
bb3=4;
bb4=8;
bb5=2;
bb6=9;
bb7=1;
vspec=Element(16,((j-1)*9)+1:((j-1)*9)+3)-Element(16,((j-1)*9)+4:((j-1)*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vector1=Element(16,((j-1)*9)+4:((j-1)*9)+6)+((Element(16,((j-1)*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
vspec=Element(8,((j-1)*9)+1:((j-1)*9)+3)-Element(8,((j-1)*9)+4:((j-1)*9)+6);
vspecn=[-vspec(1,2), vspec(1,1), 0];
vector2=Element(8,((j-1)*9)+4:((j-1)*9)+6)+((Element(8,((j-1)*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fvv=[((min-max)/2)+max, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    N=eye(6);
    N(4:6,1)=transpose(cross(vector1,[1 0 0]));
    N(4:6,2)=transpose(cross(vector1,[0 1 0]));
    N(4:6,3)=transpose(cross(vector1,[0 0 1]));
    NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
    Disp1=NewT(4,1);
    N=eye(6);
    N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
    N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
    N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
    NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
    Disp2=-NewT(4,1);
    Disp=Disp1+Disp2;
    Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
    if(Disp<Dispinc)
        Force2=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force2=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force2=((min-max)/2)+max;
        error=0;
    end
end
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fv=[0, ((min-max)/2)+max, 0];
    Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
    Fvv=[Force2, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    Disp=Twist(((6*bb1)-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
for hhh=1:HH
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, Force1, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[((min-max)/2)+max, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        N=eye(6);
        N(4:6,1)=transpose(cross(vector1,[1 0 0]));
        N(4:6,2)=transpose(cross(vector1,[0 1 0]));
        N(4:6,3)=transpose(cross(vector1,[0 0 1]));
        NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
        Disp1=NewT(4,1);
        N=eye(6);
        N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
        Disp2=-NewT(4,1);
        Disp=Disp1+Disp2;
        Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
        if(Disp<Dispinc)
            Force2=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force2=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force2=((min-max)/2)+max;
            error=0;
        end
    end
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, ((min-max)/2)+max, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[Force2, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        Disp=Twist(((6*bb1)-6)+5,1);
        Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
        if(Disp<Dispinc)
            Force1=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force1=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force1=((min-max)/2)+max;
            error=0;
        end
    end
end
Wrench=zeros(6*Constraint(21,1),1);
Fv=[0, Force1, 0];
Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
Fvv=[Force2, 0, 0];
Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS15(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS15(i,1))-6+1:(6*CS15(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS15(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS15(i,2))-6+1:(6*CS15(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS15(1:20,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase16(j,Element,CS16,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=18;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS16(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

bb1=18;
bb2=12;
bb3=13;
bb4=9;
bb5=16;
bb6=10;
bb7=17;
vspec=Element(16,((j-1)*9)+1:((j-1)*9)+3)-Element(16,((j-1)*9)+4:((j-1)*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vector1=Element(16,((j-1)*9)+4:((j-1)*9)+6)+((Element(16,((j-1)*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
vspec=Element(8,((j-1)*9)+1:((j-1)*9)+3)-Element(8,((j-1)*9)+4:((j-1)*9)+6);
vspecn=[-vspec(1,2), vspec(1,1), 0];
vector2=Element(8,((j-1)*9)+4:((j-1)*9)+6)+((Element(8,((j-1)*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fvv=[((min-max)/2)+max, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    N=eye(6);
    N(4:6,1)=transpose(cross(vector1,[1 0 0]));
    N(4:6,2)=transpose(cross(vector1,[0 1 0]));
    N(4:6,3)=transpose(cross(vector1,[0 0 1]));
    NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
    Disp1=NewT(4,1);
    N=eye(6);
    N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
    N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
    N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
    NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
    Disp2=-NewT(4,1);
    Disp=Disp1+Disp2;
    Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
    if(Disp<Dispinc)
        Force2=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force2=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force2=((min-max)/2)+max;
        error=0;
    end
end
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fv=[0, ((min-max)/2)+max, 0];
    Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
    Fvv=[Force2, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    Disp=Twist(((6*bb1)-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
for hhh=1:HH
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, Force1, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[((min-max)/2)+max, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        N=eye(6);
        N(4:6,1)=transpose(cross(vector1,[1 0 0]));
        N(4:6,2)=transpose(cross(vector1,[0 1 0]));
        N(4:6,3)=transpose(cross(vector1,[0 0 1]));
        NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
        Disp1=NewT(4,1);
        N=eye(6);
        N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
        Disp2=-NewT(4,1);
        Disp=Disp1+Disp2;
        Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
        if(Disp<Dispinc)
            Force2=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force2=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force2=((min-max)/2)+max;
            error=0;
        end
    end
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, ((min-max)/2)+max, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[Force2, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        Disp=Twist(((6*bb1)-6)+5,1);
        Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
        if(Disp<Dispinc)
            Force1=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force1=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force1=((min-max)/2)+max;
            error=0;
        end
    end
end
Wrench=zeros(6*Constraint(21,1),1);
Fv=[0, Force1, 0];
Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
Fvv=[Force2, 0, 0];
Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS16(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS16(i,1))-6+1:(6*CS16(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS16(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS16(i,2))-6+1:(6*CS16(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS16(1:20,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase17(j,Element,CS17,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=19;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS17(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

bb1=19;
bb2=13;
bb3=14;
bb4=3;
bb5=1;
bb6=4;
bb7=2;
vv1=Element(16,((j-1)*9)+1:((j-1)*9)+3)-Element(14,((j-1)*9)+1:((j-1)*9)+3);
vector1=Element(20,((j-1)*9)+1:((j-1)*9)+3)+[(Element(20,((j-1)*9)+7)/2)*sin(theta1), (Element(20,((j-1)*9)+7)/2)*cos(theta1), 0]+(((Element(20,((j-1)*9)+7)*1.5)-(Element(20,((j-1)*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
vv2=Element(8,((j-1)*9)+1:((j-1)*9)+3)-Element(7,((j-1)*9)+1:((j-1)*9)+3);
vector2=Element(2,((j-1)*9)+1:((j-1)*9)+3)+[(Element(2,((j-1)*9)+7)/2)*sin(theta1), -(Element(2,((j-1)*9)+7)/2)*cos(theta1), 0]+(((Element(2,((j-1)*9)+7)*1.5)-(Element(2,((j-1)*9)+7)*sin(theta1)))*vv2/sqrt(dot(vv2,vv2)));

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fvv=[((min-max)/2)+max, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    N=eye(6);
    N(4:6,1)=transpose(cross(vector1,[1 0 0]));
    N(4:6,2)=transpose(cross(vector1,[0 1 0]));
    N(4:6,3)=transpose(cross(vector1,[0 0 1]));
    NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
    Disp1=NewT(4,1);
    N=eye(6);
    N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
    N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
    N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
    NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
    Disp2=-NewT(4,1);
    Disp=Disp1+Disp2;
    Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
    if(Disp<Dispinc)
        Force2=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force2=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force2=((min-max)/2)+max;
        error=0;
    end
end
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fv=[0, ((min-max)/2)+max, 0];
    Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
    Fvv=[Force2, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    Disp=Twist(((6*bb1)-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
for hhh=1:HH
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, Force1, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[((min-max)/2)+max, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        N=eye(6);
        N(4:6,1)=transpose(cross(vector1,[1 0 0]));
        N(4:6,2)=transpose(cross(vector1,[0 1 0]));
        N(4:6,3)=transpose(cross(vector1,[0 0 1]));
        NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
        Disp1=NewT(4,1);
        N=eye(6);
        N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
        Disp2=-NewT(4,1);
        Disp=Disp1+Disp2;
        Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
        if(Disp<Dispinc)
            Force2=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force2=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force2=((min-max)/2)+max;
            error=0;
        end
    end
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, ((min-max)/2)+max, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[Force2, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        Disp=Twist(((6*bb1)-6)+5,1);
        Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
        if(Disp<Dispinc)
            Force1=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force1=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force1=((min-max)/2)+max;
            error=0;
        end
    end
end
Wrench=zeros(6*Constraint(21,1),1);
Fv=[0, Force1, 0];
Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
Fvv=[Force2, 0, 0];
Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS17(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS17(i,1))-6+1:(6*CS17(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS17(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS17(i,2))-6+1:(6*CS17(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS17(1:20,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase18(j,Element,CS18,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=18;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS18(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

bb1=18;
bb2=12;
bb3=13;
bb4=3;
bb5=1;
bb6=4;
bb7=2;
vv1=Element(16,((j-1)*9)+1:((j-1)*9)+3)-Element(14,((j-1)*9)+1:((j-1)*9)+3);
vector1=Element(20,((j-1)*9)+1:((j-1)*9)+3)+[(Element(20,((j-1)*9)+7)/2)*sin(theta1), (Element(20,((j-1)*9)+7)/2)*cos(theta1), 0]+(((Element(20,((j-1)*9)+7)*1.5)-(Element(20,((j-1)*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
vv2=Element(8,((j-1)*9)+1:((j-1)*9)+3)-Element(7,((j-1)*9)+1:((j-1)*9)+3);
vector2=Element(2,((j-1)*9)+1:((j-1)*9)+3)+[(Element(2,((j-1)*9)+7)/2)*sin(theta1), -(Element(2,((j-1)*9)+7)/2)*cos(theta1), 0]+(((Element(2,((j-1)*9)+7)*1.5)-(Element(2,((j-1)*9)+7)*sin(theta1)))*vv2/sqrt(dot(vv2,vv2)));

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fvv=[((min-max)/2)+max, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    N=eye(6);
    N(4:6,1)=transpose(cross(vector1,[1 0 0]));
    N(4:6,2)=transpose(cross(vector1,[0 1 0]));
    N(4:6,3)=transpose(cross(vector1,[0 0 1]));
    NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
    Disp1=NewT(4,1);
    N=eye(6);
    N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
    N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
    N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
    NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
    Disp2=-NewT(4,1);
    Disp=Disp1+Disp2;
    Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
    if(Disp<Dispinc)
        Force2=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force2=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force2=((min-max)/2)+max;
        error=0;
    end
end
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fv=[0, ((min-max)/2)+max, 0];
    Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
    Fvv=[Force2, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    Disp=Twist(((6*bb1)-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
for hhh=1:HH
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, Force1, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[((min-max)/2)+max, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        N=eye(6);
        N(4:6,1)=transpose(cross(vector1,[1 0 0]));
        N(4:6,2)=transpose(cross(vector1,[0 1 0]));
        N(4:6,3)=transpose(cross(vector1,[0 0 1]));
        NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
        Disp1=NewT(4,1);
        N=eye(6);
        N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
        Disp2=-NewT(4,1);
        Disp=Disp1+Disp2;
        Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
        if(Disp<Dispinc)
            Force2=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force2=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force2=((min-max)/2)+max;
            error=0;
        end
    end
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, ((min-max)/2)+max, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[Force2, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        Disp=Twist(((6*bb1)-6)+5,1);
        Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
        if(Disp<Dispinc)
            Force1=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force1=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force1=((min-max)/2)+max;
            error=0;
        end
    end
end
Wrench=zeros(6*Constraint(21,1),1);
Fv=[0, Force1, 0];
Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
Fvv=[Force2, 0, 0];
Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS18(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS18(i,1))-6+1:(6*CS18(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS18(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS18(i,2))-6+1:(6*CS18(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS18(1:20,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase19(j,Element,CS19,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=10;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS19(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force1; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS19(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS19(i,1))-6+1:(6*CS19(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS19(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS19(i,2))-6+1:(6*CS19(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS19(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase20(j,Element,CS20,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=9;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS20(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force1; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS20(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS20(i,1))-6+1:(6*CS20(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS20(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS20(i,2))-6+1:(6*CS20(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS20(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase21(j,Element,CS21,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=6;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS21(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force1; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS21(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS21(i,1))-6+1:(6*CS21(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS21(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS21(i,2))-6+1:(6*CS21(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS21(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase22(j,Element,CS22,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=5;
Constraint(21,3)=20;
for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS22(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; ((min-max)/2)+max; 0; 0; 0; 0];
    Twist=K\Wrench;
    Disp=Twist(((6*Constraint(21,1))-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
Wrench(((6*Constraint(21,1))-6)+1:(6*Constraint(21,1)),1)=[0; Force1; 0; 0; 0; 0];

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS22(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS22(i,1))-6+1:(6*CS22(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS22(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS22(i,2))-6+1:(6*CS22(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS22(:,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;
function [Element,strainPoisson,g1,g3,g4,g5,g6]=CCase23(j,Element,CS23,strainload,strainPoisson,la,theta1,theta2,W,dim,Vertic,Horiz,E,G,n,HH)
% Find the stiffness matrix for the cell in its current configuration
Constraint=zeros(21,16);
Constraint(21,1)=15;
Constraint(21,3)=20;

for i=1:20
    n3(i,1:3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)-Element(i,((j-1)*9)+4:((j-1)*9)+6);
    Constraint(i,1:14)=[Element(i,((j-1)*9)+1:((j-1)*9)+3), n3(i,1:3), [n3(i,2), -n3(i,1), 0], sqrt(dot(n3(i,1:3),n3(i,1:3))), W, Element(i,((j-1)*9)+7), E, G];
    Constraint(i,15:16)=CS23(i,1:2);
end
[K]=EulerStiffnessMatrixConnect(Constraint);

bb1=15;
bb2=6;
bb3=5;
bb4=12;
bb5=2;
bb6=11;
bb7=1;
vv1=Element(16,((j-1)*9)+1:((j-1)*9)+3)-Element(14,((j-1)*9)+1:((j-1)*9)+3);
vector1=Element(20,((j-1)*9)+1:((j-1)*9)+3)+[(Element(20,((j-1)*9)+7)/2)*sin(theta1), (Element(20,((j-1)*9)+7)/2)*cos(theta1), 0]+(((Element(20,((j-1)*9)+7)*1.5)-(Element(20,((j-1)*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
vv2=Element(8,((j-1)*9)+1:((j-1)*9)+3)-Element(7,((j-1)*9)+1:((j-1)*9)+3);
vector2=Element(2,((j-1)*9)+1:((j-1)*9)+3)+[(Element(2,((j-1)*9)+7)/2)*sin(theta1), -(Element(2,((j-1)*9)+7)/2)*cos(theta1), 0]+(((Element(2,((j-1)*9)+7)*1.5)-(Element(2,((j-1)*9)+7)*sin(theta1)))*vv2/sqrt(dot(vv2,vv2)));

% Determine the force required to incrementally load the cell properly
Force=-E*la*(dim^2);
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fvv=[((min-max)/2)+max, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    N=eye(6);
    N(4:6,1)=transpose(cross(vector1,[1 0 0]));
    N(4:6,2)=transpose(cross(vector1,[0 1 0]));
    N(4:6,3)=transpose(cross(vector1,[0 0 1]));
    NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
    Disp1=NewT(4,1);
    N=eye(6);
    N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
    N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
    N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
    NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
    Disp2=-NewT(4,1);
    Disp=Disp1+Disp2;
    Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
    if(Disp<Dispinc)
        Force2=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force2=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force2=((min-max)/2)+max;
        error=0;
    end
end
max=0;
min=Force;
Wrench=zeros(6*Constraint(21,1),1);
error=1e10;
while(error>1e-10)
    Fv=[0, ((min-max)/2)+max, 0];
    Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
    Fvv=[Force2, 0, 0];
    Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
    Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
    Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
    Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
    Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
    Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
    Twist=K\Wrench;
    Disp=Twist(((6*bb1)-6)+5,1);
    Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
    if(Disp<Dispinc)
        Force1=((min-max)/2)+max;
        min=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    elseif(Disp>Dispinc)
        Force1=((min-max)/2)+max;
        max=((min-max)/2)+max;
        error=abs(Disp-Dispinc);
    else
        Force1=((min-max)/2)+max;
        error=0;
    end
end
for hhh=1:HH
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, Force1, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[((min-max)/2)+max, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        N=eye(6);
        N(4:6,1)=transpose(cross(vector1,[1 0 0]));
        N(4:6,2)=transpose(cross(vector1,[0 1 0]));
        N(4:6,3)=transpose(cross(vector1,[0 0 1]));
        NewT=N\Twist((6*bb4)-6+1:(6*bb4),1);
        Disp1=NewT(4,1);
        N=eye(6);
        N(4:6,1)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*bb2)-6+1:(6*bb2),1);
        Disp2=-NewT(4,1);
        Disp=Disp1+Disp2;
        Dispinc=Element(12,((j-1)*9)+1)-vector1(1,1);
        if(Disp<Dispinc)
            Force2=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force2=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force2=((min-max)/2)+max;
            error=0;
        end
    end
    max=0;
    min=Force;
    Wrench=zeros(6*Constraint(21,1),1);
    error=1e10;
    while(error>1e-10)
        Fv=[0, ((min-max)/2)+max, 0];
        Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
        Fvv=[Force2, 0, 0];
        Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
        Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
        Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
        Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
        Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
        Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);
        Twist=K\Wrench;
        Disp=Twist(((6*bb1)-6)+5,1);
        Dispinc=(strainload(1,j+1)-strainload(1,j))*(2*Vertic);
        if(Disp<Dispinc)
            Force1=((min-max)/2)+max;
            min=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        elseif(Disp>Dispinc)
            Force1=((min-max)/2)+max;
            max=((min-max)/2)+max;
            error=abs(Disp-Dispinc);
        else
            Force1=((min-max)/2)+max;
            error=0;
        end
    end
end
Wrench=zeros(6*Constraint(21,1),1);
Fv=[0, Force1, 0];
Wrench(((6*bb1)-6)+1:(6*bb1),1)=transpose([Fv, 0, 0, 0]);
Fvv=[Force2, 0, 0];
Wrench(((6*bb2)-6)+1:(6*bb2),1)=transpose([-Fvv, cross(Element(12,((j-1)*9)+1:((j-1)*9)+3),-Fvv)]);
Wrench(((6*bb3)-6)+1:(6*bb3),1)=transpose([Fvv, cross(Element(9,((j-1)*9)+1:((j-1)*9)+3),Fvv)]);
Wrench(((6*bb4)-6)+1:(6*bb4),1)=transpose([Fvv/2, cross(vector1,Fvv/2)]);
Wrench(((6*bb5)-6)+1:(6*bb5),1)=transpose([Fvv/2, cross(vector2,Fvv/2)]);
Wrench(((6*bb6)-6)+1:(6*bb6),1)=transpose([-Fvv/2, cross([-vector1(1,1), vector1(1,2), 0],-Fvv/2)]);
Wrench(((6*bb7)-6)+1:(6*bb7),1)=transpose([-Fvv/2, cross([-vector2(1,1), vector2(1,2), 0],-Fvv/2)]);

% Update the Element matrix now that the correct load has been identified
Twist=K\Wrench;
N=eye(6);
for i=1:20
    if(CS23(i,1)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+1:((j-1)*9)+3),[0 0 1]));
        NewT=N\Twist((6*CS23(i,1))-6+1:(6*CS23(i,1)),1);
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+1:(j*9)+3)=Element(i,((j-1)*9)+1:((j-1)*9)+3);
    end
    if(CS23(i,2)~=0)
        N(4:6,1)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[1 0 0]));
        N(4:6,2)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 1 0]));
        N(4:6,3)=transpose(cross(Element(i,((j-1)*9)+4:((j-1)*9)+6),[0 0 1]));
        NewT=N\Twist((6*CS23(i,2))-6+1:(6*CS23(i,2)),1);
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6)+transpose(NewT(4:6,1));
    else
        Element(i,(j*9)+4:(j*9)+6)=Element(i,((j-1)*9)+4:((j-1)*9)+6);
    end
end
Element(1:20,(j*9)+7)=Element(1:20,((j-1)*9)+7);
Element(1:20,(j*9)+8:(j*9)+9)=CS23(1:20,1:2);

g1=Element(11,(j*9)+4)-Element(10,(j*9)+4);
g3=Element(15,(j*9)+5)-(Element(9,(j*9)+2)+(Element(9,(j*9)+7)/2));
vv1=Element(16,(j*9)+1:(j*9)+3)-Element(14,(j*9)+1:(j*9)+3);
vv2=Element(20,(j*9)+1:(j*9)+3)+[(Element(20,(j*9)+7)/2)*sin(theta1), (Element(20,(j*9)+7)/2)*cos(theta1), 0]+(((Element(20,(j*9)+7)*1.5)-(Element(20,(j*9)+7)*sin(theta1)))*vv1/sqrt(dot(vv1,vv1)));
g4=Element(12,(j*9)+1)-vv2(1,1);
g5=Element(19,(j*9)+5)+((Element(19,(j*9)+7)/2)/cos(theta1))-vv2(1,2);
vspec=Element(16,(j*9)+1:(j*9)+3)-Element(16,(j*9)+4:(j*9)+6);
vspecn=[vspec(1,2), -vspec(1,1), 0];
vv3=Element(16,(j*9)+4:(j*9)+6)+((Element(16,(j*9)+7)/2)*(vspecn/sqrt(dot(vspecn,vspecn))));
g6=Element(12,(j*9)+1)-vv3(1,1);

% Find the Poisson's ratio strain
strainPoisson(1,j+1)=(Element(12,(j*9)+1)-Element(12,1))/Horiz;



function nn_Callback(hObject, eventdata, handles)
% hObject    handle to nn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nn as text
%        str2double(get(hObject,'String')) returns contents of nn as a double


% --- Executes during object creation, after setting all properties.
function nn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function S5_Callback(hObject, eventdata, handles)
% hObject    handle to S5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of S5 as text
%        str2double(get(hObject,'String')) returns contents of S5 as a double


% --- Executes during object creation, after setting all properties.
function S5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to S5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function S6_Callback(hObject, eventdata, handles)
% hObject    handle to S6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of S6 as text
%        str2double(get(hObject,'String')) returns contents of S6 as a double


% --- Executes during object creation, after setting all properties.
function S6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to S6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
