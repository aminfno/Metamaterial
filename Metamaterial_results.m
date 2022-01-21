%    Copyright 2022 Flexible Research Group @ UCLA
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%      http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
function varargout = Metamaterial_results(varargin)
% METAMATERIAL_RESULTS MATLAB code for Metamaterial_results.fig
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Metamaterial_results_OpeningFcn, ...
    'gui_OutputFcn',  @Metamaterial_results_OutputFcn, ...
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

function Metamaterial_results_OpeningFcn(hObject, ~, handles, varargin)
load data.mat  logo1 plot_label%Load images data

%% GUI images
axes(handles.LOGO) % Load the logo and show it on the GUI
imshow(logo1);
grid off
set(handles.LOGO, 'Box', 'off')
handles.LOGO.XRuler.Axle.LineStyle = 'none';
handles.LOGO.YRuler.Axle.LineStyle = 'none';

axes(handles.plot_label) % Load the plot label image and show it on the GUI
imshow(plot_label);
grid off
set(handles.plot_label, 'Box', 'off')
handles.plot_label.XRuler.Axle.LineStyle = 'none';
handles.plot_label.YRuler.Axle.LineStyle = 'none';

axes(handles.animation) %  show the cell animation on the GUI
grid off
set(handles.animation, 'Box', 'off')
handles.animation.XRuler.Axle.LineStyle = 'none';
handles.animation.YRuler.Axle.LineStyle = 'none';
%axesHandle = findobj('Tag', 'animation');

global LOGIC
LOGIC= 1;
PlayAnim(hObject, handles)

handles.output = hObject;
guidata(hObject, handles);




% --- Outputs from this function are returned to the command line.
function varargout = Metamaterial_results_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
global LOGIC
if LOGIC==0
    close
    Metamaterial
end



% --- Executes on button press in saveGIF.
function saveGIF_Callback(hObject, eventdata, handles)
%msg= msgbox(' Saving the results as a GIF file ...','Please Wait');
GenerateGIF()


% --- Executes on button press in GoBack.
function GoBack_Callback(hObject, eventdata, handles)
global LOGIC
LOGIC= 0;
counter =1;
save data counter -append

function PlayAnim(hObject,handles)
global LOGIC
%%create animation
load data.mat images strainload time totalstrainPoisson
%handles.images=images;
time= time/60;
totalstrainPoisson=totalstrainPoisson * 100;
strainload= strainload*100;

%% subplot 2/ Results
set(gcf,'color','w')

yl=([min(strainload)-5 max(strainload)+5]);

hold off;
axes(handles.plot) %  show the animation of the plot on the GUI
a = plot(time,totalstrainPoisson,'b','LineWidth',2.5);
xlabel('Time (min)')
ylabel('Strain (%)')
hold on;
b = plot(time,strainload,'Color',[.5 .5 .5],'LineWidth',2.5);
hold on;
p = plot(time(1)*[1 1],yl,'k');
yyaxis right
pr=-1*totalstrainPoisson./strainload;
ipr = filloutliers(pr,'linear');
idx=~(abs(ipr)<.0000001);
L=max(abs(pr));
ylim([-1.1*L,1.1*L])
plot(time(idx),ipr(idx),'Color',[0.4940, 0.1840, 0.5560],'LineWidth',2.5,'LineStyle','--','DisplayName','Poisson''s ratio');
ylabel('Poisson''s ratio (-)')
hold off;
set(gca,'FontSize',20);


while LOGIC==1
    k=1;
    while k <= (length(time)-1)
        %% subplot 1 Unit cell shapes
        if LOGIC==0 || LOGIC==2
            return
        else
            p.XData = time(k)*[1 1];
            imshow(images{k},'Parent',handles.animation)
            pause(.05)
            k=k+1;
        end
    end
end


function GenerateGIF()
global LOGIC
LOGIC=2;
load data.mat images strainload time totalstrainPoisson logo1 plot_label

Filename = sprintf('%s.gif', datestr(now,'mm-dd-yyyy HH-MM'));
[file,path] = uiputfile(Filename);
fullname=fullfile(path,file);

%Checks if in dialog box cancel is pressed, stops generating the gif
if file==0
    LOGIC=1;
    return
end
%%Waitbar
Progressbar = waitbar(0,'','Name',sprintf('    Saving   "%s"     ',file),...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)','WindowStyle','modal');
setappdata(Progressbar,'canceling',0);
%%%

time= time/60;
totalstrainPoisson=totalstrainPoisson * 100;
strainload= strainload*100;

%% subplot 2: Results
f = figure('visible','off');
set(gcf, 'Position', get(0, 'Screensize')/1.1, 'color','w');
imshow(logo1,'Parent',subplot(18,2,[2,4]));
imshow(plot_label,'Parent',subplot(18,2,32:2:36));
yl= ([min(strainload)-5 max(strainload)+5]);
subplot(18,2,8:2:26);
set(gcf,'color','w');
hold off;
a = plot(time,totalstrainPoisson,'b','LineWidth',2);
xlabel('Time (min)')
ylabel('Strain (%)')
hold on;
b = plot(time,strainload,'Color',[.5 .5 .5],'LineWidth',2);
hold on;
p = plot(time(1)*[1 1],yl,'k');

yyaxis right
pr=-1*totalstrainPoisson./strainload;
ipr = filloutliers(pr,'linear');
idx=~(abs(ipr)<.0000001);
L=max(abs(pr));
ylim([-1.1*L,1.1*L])
plot(time(idx),ipr(idx),'Color',[0.4940, 0.1840, 0.5560],'LineWidth',2.5,'LineStyle','--','DisplayName','Poisson''s ratio');
% ylim([-1,1])
ylabel('Poisson''s ratio (-)')
hold off;
set(gca,'FontSize',20);

hold off;
p.XData = time(2)*[1 1];
for k = 1:(length(time)-1)
    set(0,'CurrentFigure',f);
    p.XData = time(k)*[1 1];
    %% subplot 1: Unit cell shapes
    imshow(images{k},'Parent',subplot(18,2,1:2:35));
    % Write to the GIF File
  


    frame = getframe(f);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,fullname,'gif', 'Loopcount',inf,'DelayTime',10/(length(time)-1));
    else
        imwrite(imind,cm,fullname,'gif','WriteMode','append','DelayTime',10/(length(time)-1));
    end
    
    %%Updating Waitbar
    if getappdata(Progressbar,'canceling')
        delete(Progressbar)
        LOGIC=1;
        return
    end
    waitbar(k/(length(time)-1),Progressbar,sprintf('%.0f%s',k*100/(length(time)-1),'% Complete'))
end
LOGIC=1;
delete(Progressbar)
