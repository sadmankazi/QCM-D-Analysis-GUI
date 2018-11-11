function varargout = QCMDAnalyze(varargin)

% Copyright (C) 2018 Kazi Sadman (Shull Research Group, Northwestern Uni.)
%
% This is Version 1 of the GUI "QCMDAnalyze"
% Source code: https://github.com/sadmankazi/QCM-D-Analysis-GUI
% See READ ME for referencing and citation recommendations
% Bug reporting should be done through GitHub
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>

% QCMDANALYZE MATLAB code for QCMDAnalyze.fig
%      QCMDANALYZE, by itself, creates a new QCMDANALYZE or raises the existing
%      singleton*.
%
%      H = QCMDANALYZE returns the handle to a new QCMDANALYZE or the handle to
%      the existing singleton*.
%
%      QCMDANALYZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QCMDANALYZE.M with the given input arguments.
%
%      QCMDANALYZE('Property','Value',...) creates a new QCMDANALYZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before QCMDAnalyze_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to QCMDAnalyze_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QCMDAnalyze

% Last Modified by GUIDE v2.5 10-Nov-2018 22:18:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @QCMDAnalyze_OpeningFcn, ...
    'gui_OutputFcn',  @QCMDAnalyze_OutputFcn, ...
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

% End initialization code - DO NOT EDIT

% --- Executes just before QCMDAnalyze is made visible.
function QCMDAnalyze_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QCMDAnalyze (see VARARGIN)

% set(gcf,'Units','Pixels','Position',get(0,'ScreenSize')) % Make gui full size
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Choose default command line output for QCMDAnalyze
handles.output = hObject;

% set(0,'DefaultAxesColorOrder',[0 0 1; 1 0 0; 0 0.5 0])
set(0,'defaulttextfontsize',18)
set(0,'defaultlinemarkersize',10);
set(0,'defaultaxeslinewidth',1.75);
set(0,'defaultpatchlinewidth',1.75);
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',2)

set(hObject,'toolbar','figure');

guidata(hObject,handles)
handles.output = hObject;

% --- Outputs from this function are returned to the command line.
function varargout = QCMDAnalyze_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadQCM.
function loadQCM_Callback(hObject, eventdata, handles)

% In this function we raw import data, clean it, and plot it before solving

[FileName,PathName] = uigetfile('*.xlsx','File Selector');
handles.filename = FileName;

% If someone cancels out of the uigetfile dialog, filename and pathname will
% both be 0. This checks if that is the case.
if ~FileName
    set(handles.statusupdate, 'String', 'No QCM data loaded!','Foregroundcolor','red');
    return
end

handles.pathname=PathName; %Save this pathname to handle structure so we can open this later for current data

fullpathname = fullfile(PathName,FileName);
[pathstr,name,ext] = fileparts(fullpathname);
handles.filename=name; % only save the name of the file

set(handles.textQCM, 'String', FileName);
guidata(hObject, handles)

data = xlsread(fullpathname);

%Delete any missing data in the matrix:
[row, col] = find(isnan(data));
data(row, :) = [];

%Import data into the handle structure:
handles.t = data(:,1)./60; %Convert seconds to min
handles.delf{1} = data(:,2);
handles.delg{1} = 0.5*1*5*data(:,3);

handles.delf{3} = data(:,4);
handles.delg{3} = 0.5*3*5*data(:,5);

handles.delf{5} = data(:,6);
handles.delg{5} = 0.5*5*5*data(:,7);
%handles.temp = data(:,8);  Import new variables here if (temperature for
%example), put these data into the 8th, 9th columns...

set(handles.endindex,'String',num2str(length(handles.t)));
endidx = str2num(get(handles.endindex,'String'));

plot(handles.axes4,  handles.t, handles.delf{3}./3,'o',...
    handles.t, handles.delf{5}./5,'x');
plot(handles.axes6, handles.t, handles.delg{3},'o',...
    handles.t, handles.delg{5},'x');

ylabel(handles.axes4,'\Deltaf_n/n (Hz)','fontweight','bold');
xlabel(handles.axes4,'t (min)','fontweight','bold');
xlim(handles.axes4,[0 handles.t(endidx)]);
legend(handles.axes4,'n=3','n=5','location','best');

ylabel(handles.axes6,'\Delta\Gamma_n (Hz)','fontweight','bold');
xlabel(handles.axes6,'t (min)','fontweight','bold');
xlim(handles.axes6,[0 handles.t(endidx)]);
linkaxes([handles.axes4,handles.axes6],'x');

set(handles.statusupdate, 'String', 'QCM-D Data loaded!','Foregroundcolor',[0 0.5 0]);

guidata(hObject,handles)

% --- Executes on button press in plotqcmdata.
function plotqcmdata_Callback(hObject, eventdata, handles)

if ~isfield(handles,'filename')
    set(handles.statusupdate, 'String', 'Load QCM data first!','Foregroundcolor','red');
    return
end

%Check if shifts are liekly in the bulk regime or not:
if -handles.delf{1} < handles.delg{1}
    set(handles.statusupdate, 'String', 'Shifts not in thin-film limit! Try Bulk.','Foregroundcolor','red');
    return
end

%Check if last data index entered is larger than number of data points:
if str2num(get(handles.endindex,'String')) > length(handles.t)
    set(handles.statusupdate, 'String', 'Data range cannot exceed total no. of points!','Foregroundcolor','red');
    set(handles.endindex,'String',num2str(length(handles.t)));
    return
end

handles.edissratio = [];
handles.eharmratio = [];
handles.d1out = [];
handles.drhoout = [];
handles.phiout = [];

f1=5e6;        % Fundamental frequency
zq=8.84e6;     % Modulus of quartz, kg/m^2-s
dgliq3=1300;   % dissipation at n=3 for water

sauerbrey=@(n,drho) 2*n*f1^2*drho/zq;
delfstarliq=@(n,etarho) (f1^1.5/zq)*(n*etarho/pi)^0.5*(1i-1);

ndrho=3; % harmonic used for thickness calculation

% nh{1}=[1 3 3];
nh{1}=[3 5 3];
nh{2}=[3 5 5];
% nh{4}=[5 7 5];  %only if code is modified to import 7th harmonic above

i=length(nh);

if get(handles.onelayer,'value')==0
    handles.etarhoexpt=pi*dgliq3^2*zq^2/(3*f1^3);
    Rliq=@(n,drho) delfstarliq(n,handles.etarhoexpt)./sauerbrey(n,drho);  % defined in Elizabth's 2-layer paper
elseif get(handles.onelayer,'value')==1
    handles.etarhoexpt = 0;
    Rliq=@(n,drho) 0+0i;
end

d=@(n,d1,phi) d1.*n.^(1-phi./180);    % d/lambda
Dn=@(n,d1,phi)   2.*pi.*d(n,d1,phi).*(1-1i.*tand(phi./2));    % defined in Elizabth's 2015 Viscoelastic paper
grho=@(n,d1,drho,phi) n.^2.*f1^2.*(cosd(phi./2)).^2.*(1./d(n,d1,phi)).^2.*drho.^2;
delfstardn=@(Dn,Rliq)  -(Dn.^-2+Rliq.^2)./((cot(Dn)./Dn)+Rliq);
delfstar2layer=@(n,d1,phi,drho) delfstardn(Dn(n,d1,phi),Rliq(n,drho));
drhocalc=@(n,df,d1,phi,drho) df.*(zq./(2*n*f1^2))./real(delfstar2layer(n,d1,phi,drho));
rhodelta=@(n,grho,phi) ((grho.^0.5)./(2*pi.*n.*f1.*sind(phi./2)));

soln=[0.05, 45, 0.001];  % [not sure d1/lam/ phase angle/ Drho] These are the initial guesses for d1/lam, phi and drho.
lb=[0.01, 0, 2e-4];      % lower bounds on final solution
ub=[0.25, 90, 0.01];     % upper bounds on final solution

inputsoln= soln;

start = str2num(get(handles.startindex,'String'));
step = str2num(get(handles.stepindex,'String'));
endidx = str2num(get(handles.endindex,'String'));

options = optimset('display','off','TolFun',10e-6); % Set the solver tolerance

for i = 1:i
    
    set(handles.statusupdate, 'String', [num2str(round(i/3*100)) ' % complete.'],'Foregroundcolor',[0 0.5 0]);
    drawnow
    
    for k=start:step:endidx
        
        % a, b and c correspond to the n1:n2,n3 calculation
        dfa=handles.delf{nh{i}(1)}(k); dfb=handles.delf{nh{i}(2)}(k); dfc=handles.delf{nh{i}(3)}(k); dgc=handles.delg{nh{i}(3)}(k);
        
        handles.eharmratio(k)=nh{i}(1)*dfb/(nh{i}(2)*dfa);   %Experimental harmonic ratio
        handles.edissratio(k)=-dgc/dfc;      %Experimantal dissipation ratio
        
        % generalzed to treat n1:n2,n3 calculation
        fdissratio=@(d1,phi,drho) -imag(delfstar2layer(nh{i}(3),d1,phi,drho))/(real(delfstar2layer(nh{i}(3),d1,phi,drho)));
        fharmratio=@(d1,phi,drho) real(delfstar2layer(nh{i}(2),d1,phi,drho))/(real(delfstar2layer(nh{i}(1),d1,phi,drho)));
        
        if get(handles.twolayer,'value')==1
            f3tosolve=@(x) [fharmratio(x(1),x(2),x(3))-handles.eharmratio(k);...
                fdissratio(x(1),x(2),x(3))-handles.edissratio(k);...
                100*(drhocalc(ndrho, handles.delf{ndrho}(k),x(1),x(2),x(3))-x(3))];  %multiply by 100 to get enough accuracy
            inputsoln= soln;                               % Dynamic guesses for solving, i.e., previous solution is new initial guess
        elseif get(handles.onelayer,'value')==1
            f3tosolve=@(x) [fharmratio(x(1),x(2),0)-handles.eharmratio(k);...
                fdissratio(x(1),x(2),0)-handles.edissratio(k)];
            soln=soln(1:2);  % [not sure d1/lam/ phase angle] These are the initial guesses for d1/lam, phi and drho.
            inputsoln= soln;
        end
        
        try
            if get(handles.twolayer,'value')==1
                [soln,resnorm,residual,exitflag]=lsqnonlin(f3tosolve,inputsoln,lb,ub,options);
            elseif get(handles.onelayer,'value')==1
                [soln,resnorm,residual,exitflag]=lsqnonlin(f3tosolve,inputsoln(1:2),lb(1:2),ub(1:2),options);
                soln(3) = (dfa./real(delfstardn(Dn(ndrho,soln(1),soln(2)),0)))*zq./(2*1*f1.^2);
            end
            
        catch Err
            soln = [NaN NaN NaN];
        end
        
        if ~isnan(soln(1))
            inputsoln = soln;
        end
        
        handles.d1out{i}(k) = soln(1);
        handles.phiout{i}(k) = soln(2);
        handles.drhoout{i}(k) = drhocalc(ndrho,handles.delf{ndrho}(k),soln(1),soln(2),soln(3));
        handles.exitflag{i}(k) = exitflag;
    end
end

set(handles.statusupdate, 'String', 'Solved!','Foregroundcolor',[0 0.5 0]);

% In case the we don't solve for every point, delete the ones which weren't
% solved for so that the matrix dimensions match later for plotting
if step>1
    for i = 1:1:length(nh)
        handles.d1out{i} = handles.d1out{i}(handles.d1out{i}~=0);
        handles.phiout{i} = handles.phiout{i}(handles.phiout{i}~=0);
        handles.drhoout{i} = handles.drhoout{i}(handles.drhoout{i}~=0);
    end
end

% index 2 here means we are using th 3:5,5 solutions
handles.d3out = d(3, handles.d1out{2}, handles.phiout{2});
handles.d5out = d(5, handles.d1out{2}, handles.phiout{2});

handles.grho3out{1} = grho(3, handles.d1out{1}, handles.drhoout{1}, handles.phiout{1});
handles.grho3out{2} = grho(3, handles.d1out{2}, handles.drhoout{2}, handles.phiout{2});
handles.grho5out = grho(5, handles.d1out{2}, handles.drhoout{2}, handles.phiout{2});

% handles.sauerbreycorrection1 = real(delfstar2layer(1, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}));
handles.sauerbreycorrection3 = real(delfstar2layer(3, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}));
handles.sauerbreycorrection5 = real(delfstar2layer(5, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}));

% handles.rhodel1out = rhodelta(1, handles.grho1out, handles.phiout{2});
handles.rhodel3out = rhodelta(3, handles.grho3out{2}, handles.phiout{2});
handles.rhodel5out = rhodelta(5, handles.grho5out, handles.phiout{2});

% Calculate predicted shifts based on phi, d/lambda, drho that was just
% solved for:

f3pred355 = (2.*3.*f1^2).*(handles.drhoout{2}.*(real(delfstar2layer(3, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));
g3pred355 = (2.*3.*f1^2).*(handles.drhoout{2}.*(imag(delfstar2layer(3, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));
f5pred355 = (2.*5.*f1^2).*(handles.drhoout{2}.*(real(delfstar2layer(5, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));
g5pred355 = (2.*5.*f1^2).*(handles.drhoout{2}.*(imag(delfstar2layer(5, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));

f3pred353 = (2.*3.*f1^2).*(handles.drhoout{1}.*(real(delfstar2layer(3, handles.d1out{1}, handles.phiout{1}, handles.drhoout{1}))./zq));
g3pred353 = (2.*3.*f1^2).*(handles.drhoout{1}.*(imag(delfstar2layer(3, handles.d1out{1}, handles.phiout{1}, handles.drhoout{1}))./zq));
f5pred353 = (2.*5.*f1^2).*(handles.drhoout{1}.*(real(delfstar2layer(5, handles.d1out{1}, handles.phiout{1}, handles.drhoout{1}))./zq));
g5pred353 = (2.*5.*f1^2).*(handles.drhoout{1}.*(imag(delfstar2layer(5, handles.d1out{1}, handles.phiout{1}, handles.drhoout{1}))./zq));

% f1pred355 = (2.*1.*f1^2).*(handles.drhoout{2}.*(real(delfstar2layer(1, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));
% g1pred355 = (2.*1.*f1^2).*(handles.drhoout{2}.*(imag(delfstar2layer(1, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));

% Viscoelastic Plots
plot(handles.axes14, handles.t(start:step:endidx), 1e6*handles.drhoout{2}, '+', handles.t(start:step:endidx), 1e6*handles.drhoout{1}, '+');
plot(handles.axes10, handles.t(start:step:endidx), 0.001*handles.grho3out{2}, '+', handles.t(start:step:endidx), 0.001*handles.grho3out{1}, '+'); %Multiply 0.001 to convert from Pa*kg/m^3 to Pa*g/cm^3
plot(handles.axes11, handles.t(start:step:endidx), handles.phiout{2}, '+', handles.t(start:step:endidx), handles.phiout{1}, '+');
legend(handles.axes14,'3:5,5','3:5,3','location','best');

% Harmonic Plots
cla(handles.axes4)
cla(handles.axes6)

if get(handles.calc355,'value')==1
    plot(handles.axes4,...
        handles.t(start:step:endidx), handles.delf{3}(start:step:endidx)./3,'o',...
        handles.t(start:step:endidx), handles.delf{5}(start:step:endidx)./5,'x',...
        handles.t(start:step:endidx), f3pred355./3,'c.',...
        handles.t(start:step:endidx), f5pred355./5,'c.');
    legend(handles.axes4,'n=3','n=5','Pred')
    
    plot(handles.axes6,  ...
        handles.t(start:step:endidx), handles.delg{3}(start:step:endidx),'o', ...
        handles.t(start:step:endidx), handles.delg{5}(start:step:endidx),'x',...
        handles.t(start:step:endidx), g5pred355,'c.',...
        handles.t(start:step:endidx), g3pred355,'g.');
    
elseif get(handles.calc353,'value')==1
    
    plot(handles.axes4,...
        handles.t(start:step:endidx), handles.delf{3}(start:step:endidx)./3,'o',...
        handles.t(start:step:endidx), handles.delf{5}(start:step:endidx)./5,'x',...
        handles.t(start:step:endidx), f3pred353./3,'c.',...
        handles.t(start:step:endidx), f5pred353./5,'c.');
    legend(handles.axes4,'n=3','n=5','Pred')
    
    plot(handles.axes6,...
        handles.t(start:step:endidx), handles.delg{3}(start:step:endidx),'o', ...
        handles.t(start:step:endidx), handles.delg{5}(start:step:endidx),'x',...
        handles.t(start:step:endidx), g5pred353,'g.',...
        handles.t(start:step:endidx), g3pred353,'c.');
end

linkaxes([handles.axes4,handles.axes6,handles.axes14,handles.axes10,...
    handles.axes11],'x');

% Now put all the xlabels and ylabels...
ylabel(handles.axes4,'\Deltaf_n/n (Hz)','fontweight','bold');
xlabel(handles.axes4,'t (min)','fontweight','bold');
xlim(handles.axes4,[handles.t(start) handles.t(endidx)])

ylabel(handles.axes6,'\Delta\Gamma_n (Hz)','fontweight','bold');
xlabel(handles.axes6,'t (min)','fontweight','bold');
xlim(handles.axes6,[handles.t(start) handles.t(endidx)])

xlabel(handles.axes14,'t (min)','fontweight','bold');
ylabel(handles.axes14,'\Delta M_A (mg/m^2)','fontweight','bold');
xlim(handles.axes14,[handles.t(start) handles.t(endidx)])

ylabel(handles.axes10,'\rho |G*_3| (Pa-g/cm^3)','fontweight','bold');
xlabel(handles.axes10,'t (min)','fontweight','bold');
xlim(handles.axes10,[handles.t(start) handles.t(endidx)])

ylabel(handles.axes11,'\phi (deg.)','fontweight','bold');
xlabel(handles.axes11,'t (min)','fontweight','bold');
xlim(handles.axes11,[handles.t(start) handles.t(endidx)])

% Save predicted shifts into the handle structure:
handles.f3pred355=f3pred355;
handles.f5pred355=f5pred355;
handles.g3pred355=g3pred355;
handles.g5pred355=g5pred355;

handles.f3pred353=f3pred353;
handles.f5pred353=f5pred353;
handles.g3pred353=g3pred353;
handles.g5pred353=g5pred353;

zoom on
clear('start','step','endidx')
guidata(hObject, handles);


% --- Executes on button press in refresh.
function refresh_Callback(hObject, eventdata, handles)
close(gcbf)
close all
clear all
reset(0)
QCMDAnalyze
clc


% --- Executes on button press in verifysoln.
function verifysoln_Callback(hObject, eventdata, handles)

if ~isfield(handles,'d1out')
    set(handles.statusupdate, 'String', 'No solved solutions!','Foregroundcolor','red');
    return
end

start = str2num(get(handles.startindex,'String'));
step = str2num(get(handles.stepindex,'String'));
endidx = str2num(get(handles.endindex,'String'));

calcprops=figure('units','inches','Position', [2.5, 4, 12, 5]);
subplot(1,2,1);
plot(handles.t(start:step:endidx), handles.exitflag{1}(start:step:endidx), 'o')
xlabel('t (min)');
ylabel('Solver Exitflag');
title('3:5,3 Calculation')

subplot(1,2,2);
plot(handles.t(start:step:endidx), handles.exitflag{2}(start:step:endidx),'o')
xlabel('t (min)');
ylabel('Solver Exitflag');
title('3:5,5 Calculation')


% --- Executes on button press in showhandles.
function showhandles_Callback(hObject, eventdata, handles)
keyboard


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in onelayer.
function onelayer_Callback(hObject, eventdata, handles)
% hObject    handle to onelayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of onelayer


function endindex_Callback(hObject, eventdata, handles)
% hObject    handle to endindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endindex as text
%        str2double(get(hObject,'String')) returns contents of endindex as a double


% --- Executes during object creation, after setting all properties.
function endindex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function stepindex_Callback(hObject, eventdata, handles)
% hObject    handle to stepindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepindex as text
%        str2double(get(hObject,'String')) returns contents of stepindex as a double


% --- Executes during object creation, after setting all properties.
function stepindex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startindex_Callback(hObject, eventdata, handles)
% hObject    handle to startindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startindex as text
%        str2double(get(hObject,'String')) returns contents of startindex as a double


% --- Executes during object creation, after setting all properties.
function startindex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in saveplots.
function saveplots_Callback(hObject, eventdata, handles)
% hObject    handle to saveplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveplots

% --- Executes on button press in saveplots.
function twolayer_Callback(hObject, eventdata, handles)
% hObject    handle to saveplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveplots


% --- Executes on button press in calc353.
function calc353_Callback(hObject, eventdata, handles)
% hObject    handle to calc353 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calc353


% --- Executes on button press in calc355.
function calc355_Callback(hObject, eventdata, handles)
% hObject    handle to calc355 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calc355


% --- Executes on button press in exportdata.
function exportdata_Callback(hObject, eventdata, handles)

if ~isfield(handles,'d1out')
    set(handles.statusupdate, 'String', 'No solved solutions!','Foregroundcolor','red');
    return
end

s = str2num(get(handles.startindex,'String'));
st = str2num(get(handles.stepindex,'String'));
e = str2num(get(handles.endindex,'String'));

header={'time (min)','delf1 (Hz)','delg1 (Hz)','delf3','delg3','delf5','delg5',...
    'dp (mg/m2) 355','rhoG_3 355 (Pa-g/cm3)','phi 355(deg)','d/lambda_3 355','d/lambda_5 355',...
    'decay length_3 (g/m2)','decay length_5 (g/m2)','pd/(pd)_s3','delf3pred355','delg3pred355','delf5pred355','delg5pred355','Exitflag_355'};

A = [handles.t(s:st:e) handles.delf{1}(s:st:e) handles.delg{1}(s:st:e) handles.delf{3}(s:st:e) handles.delg{3}(s:st:e)...
    handles.delf{5}(s:st:e) handles.delg{5}(s:st:e) 1e6*handles.drhoout{2}' 0.001*handles.grho3out{2}'...
    handles.phiout{2}' handles.d3out' handles.d5out' 1e3*handles.rhodel3out' 1e3*handles.rhodel5out'...
    (-1./handles.sauerbreycorrection3)' handles.f3pred355' handles.g3pred355' handles.f5pred355' handles.g5pred355' handles.exitflag{2}(s:st:e)'];

fileID = fopen('exported_data.txt','w');
fprintf(fileID, '%s\t', header{:});
fclose(fileID);
dlmwrite('exported_data.txt',A,'delimiter','\t','-append','roffset',1)


% --- Executes on button press in bulkcalg.
function bulkcalc_Callback(hObject, eventdata, handles)

% Calculates pG and phi for bulk liquid measurements
%Check is any file is actually loaded:
if ~isfield(handles,'filename')
    set(handles.statusupdate, 'String', 'No QCM data loaded!','Foregroundcolor','red');
    return
end

%Check if shifts are liekly in the bulk regime or not:
if -handles.delf{1} > handles.delg{1}+100
    set(handles.statusupdate, 'String', 'Shifts not in the bulk limit!','Foregroundcolor','red');
    return
end

calcprops=figure('units','inches','Position', [2.5, 4, 12, 5]);

%Calculate phase angle based on bulk limit equations:
% phi1=2.*atan(-handles.delf{1}./handles.delg{1}).*(180/pi);
phi3=2.*atan(-handles.delf{3}./handles.delg{3}).*(180/pi);
phi5=2.*atan(-handles.delf{5}./handles.delg{5}).*(180/pi);

% Calculate Complex Modulus from phase angle
% pG1=(handles.delf{1}.*pi.*8.84e6./(5e6.*sind(phi1./2))).^2;
pG3=(handles.delf{3}.*pi.*8.84e6./(5e6.*sind(phi3./2))).^2;
pG5=(handles.delf{5}.*pi.*8.84e6./(5e6.*sind(phi5./2))).^2;

subplot(1,2,1)
plot(handles.t,pG3./1000,'o',handles.t,pG5./1000,'x');
xlabel('t (min)','fontweight','bold');
ylabel('\rho|G_n^*| (Pa-g/cm^3)','fontweight','bold');
title('(a)','fontweight','bold')

legend('n=3','n=5', 'location','best');
legend boxoff

subplot(1,2,2)
plot(handles.t,phi3,'o',handles.t,phi5,'x');
xlabel('t (min)','fontweight','bold');
ylabel('\phi (Deg.)','fontweight','bold','fontweight','bold');
title('(b)')
% ylim([0 90])

% --- Executes on button press in thicknessplots.
function thicknessplots_Callback(hObject, eventdata, handles)

%Check if solutions exist:
if ~isfield(handles,'d1out')
    set(handles.statusupdate, 'String', 'No solved solutions!','Foregroundcolor','red');
    return
end

endidx = str2num(get(handles.endindex,'String'));
start = str2num(get(handles.startindex,'String'));
step = str2num(get(handles.stepindex,'String'));

% Thickness Plots:
calcprops=figure('units','inches','Position', [2.5, 4, 17, 5]);

subplot(1,3,1)
plot(handles.t(start:step:endidx),handles.d3out,'o',handles.t(start:step:endidx),handles.d5out,'x');
ylabel('d/\lambda_n');
xlabel('t (min)');
xlim([handles.t(start) handles.t(endidx)])
ylim([0.03 0.2])
title('(a)','fontweight','bold')
legend('n=3','n=5','Location','best')
legend boxoff

subplot(1,3,2)
plot( handles.t(start:step:endidx), -1./handles.sauerbreycorrection3,'o',handles.t(start:step:endidx), -1./handles.sauerbreycorrection5,'x');
ylabel('\rhod/(\rhod)_{sn}');
xlabel('t (min)');
xlim([handles.t(start) handles.t(endidx)])
title('(b)','fontweight','bold')

subplot(1,3,3)
plot(handles.t(start:step:endidx), handles.rhodel3out*1000,'o',handles.t(start:step:endidx), handles.rhodel5out*1000,'x'); %Multiply by 1000 to convert g/m^2
ylabel('\rho\delta_n (g/m^2)');
xlabel('t (min)');
title('(c)','fontweight','bold')
xlim([handles.t(start) handles.t(endidx)])


if get(handles.saveplots,'value')==1
    calcprops.PaperPosition=[0 0 18 5];
    calcprops.PaperSize=[18 5];
    print(calcprops,'thickness plots.eps','-depsc')
    saveas(calcprops,'thickness plots.svg')
else
end


% Let's plot the frequency shifts and predicted shifts now:
calcprops=figure('units','inches','Position', [2.5, 4, 12, 12]);

if get(handles.calc355,'value')==1
    a=subplot(2,2,1);
    plot(handles.t(start:step:endidx),handles.delf{3}(start:step:endidx),'+', handles.t(start:step:endidx), handles.f3pred355,'o');
    
    b=subplot(2,2,2);
    plot(handles.t(start:step:endidx),handles.delg{3}(start:step:endidx),'+')
    hold on
    plot(handles.t(start:step:endidx), handles.g3pred355,'o','color', [0.4660    0.6740    0.1880]);
    hold off
    
    c=subplot(2,2,3);
    plot(handles.t(start:step:endidx),handles.delf{5}(start:step:endidx),'+', handles.t(start:step:endidx), handles.f5pred355,'o');
    
    d=subplot(2,2,4);
    plot(handles.t(start:step:endidx),handles.delg{5}(start:step:endidx),'+', handles.t(start:step:endidx), handles.g5pred355,'o');
    
elseif get(handles.calc353,'value')==1
    
    a=subplot(2,2,1);
    plot(handles.t(start:step:endidx),handles.delf{3}(start:step:endidx),'+', handles.t(start:step:endidx), handles.f3pred353,'o');
    
    b=subplot(2,2,2);
    plot(handles.t(start:step:endidx),handles.delg{3}(start:step:endidx),'+', handles.t(start:step:endidx), handles.g3pred353,'o')
    
    c=subplot(2,2,3);
    plot(handles.t(start:step:endidx),handles.delf{5}(start:step:endidx),'+', handles.t(start:step:endidx), handles.f5pred353,'o');
    
    d=subplot(2,2,4);
    plot(handles.t(start:step:endidx),handles.delg{5}(start:step:endidx),'+');
    hold on
    plot(handles.t(start:step:endidx),handles.g5pred353,'o','color', [0.4660    0.6740    0.1880]);
    hold off
end

ylabel(a,'\Deltaf_3 (Hz)');
xlabel(a,'t (min)');
xlim(a,[handles.t(start) handles.t(endidx)])
title(a,'(a)','fontweight','bold')
leg1=legend(a,'Experimental','Predicted','Location','best');
set(leg1,'Box','off')

ylabel(b,'\Delta\Gamma_3 (Hz)');
xlabel(b,'t (min)');
xlim(b,[handles.t(start) handles.t(endidx)])
title(b,'(b)','fontweight','bold')

ylabel(c,'\Deltaf_5 (Hz)');
xlabel(c,'t (min)');
title(c,'(c)','fontweight','bold')
xlim(c,[handles.t(start) handles.t(endidx)])

ylabel(d,'\Delta\Gamma_5 (Hz)');
xlabel(d,'t (min)');
title(d,'(d)','fontweight','bold')
xlim(d,[handles.t(start) handles.t(endidx)])

if get(handles.saveplots,'value')==1
    calcprops.PaperPosition=[0 0 12 12];
    calcprops.PaperSize=[12 12];
    print(calcprops,'shifts.eps','-depsc')
    saveas(calcprops,'shifts.svg')
else
end

% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)

if ~isfield(handles,'d1out')
    set(handles.statusupdate, 'String', 'No solved values to plot!','Foregroundcolor','red');
    return
end

calcprops=figure('units','inches','Position', [2.5, 4, 18, 5]);
% Syntax: set(gcf,?position?,[a b W H])
% (a,b) = is the lower left corner
% H = the horizontal length of the plot window
% W = the vertical width of the plot window

endidx = str2num(get(handles.endindex,'String'));
start = str2num(get(handles.startindex,'String'));
step = str2num(get(handles.stepindex,'String'));

if get(handles.calc355,'value')==1
    a = subplot(1,3,1);
    plot(a, handles.t(start:step:endidx), 1e6*handles.drhoout{2}, '+');
    
    b = subplot(1,3,2);
    semilogy(b,handles.t(start:step:endidx), 0.001*handles.grho3out{2}', '+'); %Multiply 0.001 to convert from Pa*kg/m^3 to Pa*g/cm^3
    
    c = subplot(1,3,3);
    plot(c,handles.t(start:step:endidx),handles.phiout{2}', '+');
    
elseif get(handles.calc353,'value')==1
    a = subplot(1,3,1);
    plot(a, handles.t(start:step:endidx), 1e6*handles.drhoout{1}, '+');
    
    b = subplot(1,3,2);
    semilogy(b,handles.t(start:step:endidx), 0.001*handles.grho3out{1}', '+'); %Multiply 0.001 to convert from Pa*kg/m^3 to Pa*g/cm^3
    
    c = subplot(1,3,3);
    plot(c,handles.t(start:step:endidx),handles.phiout{1}', '+');
end

xlabel(a,'t (min) ');
ylabel(a,'\rhod (mg/m^2)')
xlim(a,[handles.t(start) handles.t(endidx)])
title(a,'(a)','fontweight','bold');

ylabel(b,'Shear Modulus |G^*_3|\rho (Pa-g/cm^3)')
title(b,'(b)','fontweight','bold')
xlabel(b,'t (min)')
xlim(b,[handles.t(start) handles.t(endidx)])

xlabel(c,'t (min)');
ylabel(c,'\phi (deg.)');
title(c,'(c)','fontweight','bold');
xlim(c,[handles.t(start) handles.t(endidx)])

if get(handles.saveplots,'value')==1
    calcprops.PaperPosition=[0 0 18 5];
    calcprops.PaperSize=[18 5];
    print(calcprops,'viscoelastic plots.eps','-depsc')
    saveas(calcprops,'viscoelastic plots.svg')
else
end
