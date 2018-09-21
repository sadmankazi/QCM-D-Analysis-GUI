function varargout = shiftsimulator(varargin)
% SHIFTSIMULATOR MATLAB code for shiftsimulator.fig
%      SHIFTSIMULATOR, by itself, creates a new SHIFTSIMULATOR or raises the existing
%      singleton*.
%
%      H = SHIFTSIMULATOR returns the handle to a new SHIFTSIMULATOR or the handle to
%      the existing singleton*.
%
%      SHIFTSIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHIFTSIMULATOR.M with the given input arguments.
%
%      SHIFTSIMULATOR('Property','Value',...) creates a new SHIFTSIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before shiftsimulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to shiftsimulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help shiftsimulator

% Last Modified by GUIDE v2.5 08-Feb-2018 12:02:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @shiftsimulator_OpeningFcn, ...
                   'gui_OutputFcn',  @shiftsimulator_OutputFcn, ...
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


% --- Executes just before shiftsimulator is made visible.
function shiftsimulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to shiftsimulator (see VARARGIN)

% Choose default command line output for shiftsimulator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes shiftsimulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = shiftsimulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function rhod_Callback(hObject, eventdata, handles)
% hObject    handle to rhod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhod as text
%        str2double(get(hObject,'String')) returns contents of rhod as a double


% --- Executes during object creation, after setting all properties.
function rhod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rhog_Callback(hObject, eventdata, handles)
% hObject    handle to rhog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhog as text
%        str2double(get(hObject,'String')) returns contents of rhog as a double


% --- Executes during object creation, after setting all properties.
function rhog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phi_Callback(hObject, eventdata, handles)
% hObject    handle to phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phi as text
%        str2double(get(hObject,'String')) returns contents of phi as a double


% --- Executes during object creation, after setting all properties.
function phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delf1_Callback(hObject, eventdata, handles)
% hObject    handle to delf1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delf1 as text
%        str2double(get(hObject,'String')) returns contents of delf1 as a double


% --- Executes during object creation, after setting all properties.
function delf1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delf1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delf3_Callback(hObject, eventdata, handles)
% hObject    handle to delf3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delf3 as text
%        str2double(get(hObject,'String')) returns contents of delf3 as a double


% --- Executes during object creation, after setting all properties.
function delf3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delf3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delf5_Callback(hObject, eventdata, handles)
% hObject    handle to delf5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delf5 as text
%        str2double(get(hObject,'String')) returns contents of delf5 as a double


% --- Executes during object creation, after setting all properties.
function delf5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delf5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delg1_Callback(hObject, eventdata, handles)
% hObject    handle to delg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delg1 as text
%        str2double(get(hObject,'String')) returns contents of delg1 as a double


% --- Executes during object creation, after setting all properties.
function delg1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delg3_Callback(hObject, eventdata, handles)
% hObject    handle to delg3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delg3 as text
%        str2double(get(hObject,'String')) returns contents of delg3 as a double


% --- Executes during object creation, after setting all properties.
function delg3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delg3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delg5_Callback(hObject, eventdata, handles)
% hObject    handle to delg5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delg5 as text
%        str2double(get(hObject,'String')) returns contents of delg5 as a double


% --- Executes during object creation, after setting all properties.
function delg5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delg5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f1=5e6;        % Hz, Fundamental frequency
zq=8.84e6;     % Modulus of quartz, kg/m^2-s

sauerbrey=@(n,drho) 2*n*f1^2*drho/zq;
delfstarliq=@(n,etarho) (f1^1.5/zq)*(n*etarho/pi)^0.5*(1i-1);

dgliq =1300; %Shift in g3 only due to the liquid (assuming water)
 
if get(handles.filminliquid,'value')==1
    etarhoexpt=pi*dgliq^2*zq^2/(3*f1^3);
    Rliq=@(n,drho) delfstarliq(n,etarhoexpt)./sauerbrey(n,drho);  % defined in Elizabth's 2-layer paper
elseif get(handles.filminliquid,'value')==0
    etarhoexpt = 0;
    Rliq=@(n,drho) 0+0i;
end

d=@(n,d1,phi) d1.*n.^(1-phi./180);    % d/lambda
Dn=@(n,d1,phi)   2.*pi.*d(n,d1,phi).*(1-1i.*tand(phi./2));    % defined in Elizabth's 2015 Viscoelastic paper
delfstardn=@(Dn,Rliq)  -(Dn.^-2+Rliq.^2)./((cot(Dn)./Dn)+Rliq);
delfstar2layer=@(n,d1,phi,drho) delfstardn(Dn(n,d1,phi),Rliq(n,drho));

pd = 1e-6*str2double(get(handles.rhod,'string')); %mg/m^2
rhog3 = 1000*str2double(get(handles.rhog,'string')); %Pa-g/cm^3
phi = str2double(get(handles.phi,'string')); %deg.
d=@(n,d1,phi) d1.*n.^(1-phi./180);

n=1;
refn=3;
lambdarhon = lambdarhof(refn, n, rhog3, phi); %Lambda rho for the harmonic of interest
d1 = pd./lambdarhon;

d3 = d(3,d1,phi);
d5 = d(5,d1,phi);
d7 = d(7,d1,phi);
% d9 = d(9,d1,phi);

dev3 = 100*abs(real(delfstardn(Dn(3,d1,phi),Rliq(3,pd)))+1);
dev5 = 100*abs(real(delfstardn(Dn(5,d1,phi),Rliq(5,pd)))+1);


pdecay1 = 1e3*lambdarhof(3, 1, rhog3, phi)*cotd(phi*0.5)/(2*pi);
pdecay3 = 1e3*lambdarhof(3, 3, rhog3, phi)*cotd(phi*0.5)/(2*pi);
pdecay5 = 1e3*lambdarhof(3, 5, rhog3, phi)*cotd(phi*0.5)/(2*pi);
pdecay7 = 1e3*lambdarhof(3, 7, rhog3, phi)*cotd(phi*0.5)/(2*pi);

f1pred = (2.*1.*f1^2).*(pd*(real(delfstar2layer(1, d1, phi, pd))./zq));
f3pred = (2.*3.*f1^2).*(pd*(real(delfstar2layer(3, d1, phi, pd))./zq));
f5pred = (2.*5.*f1^2).*(pd*(real(delfstar2layer(5, d1, phi, pd))./zq));
f7pred = (2.*7.*f1^2).*(pd*(real(delfstar2layer(7, d1, phi, pd))./zq));

g1pred = (2.*1.*f1^2).*(pd*(imag(delfstar2layer(1, d1, phi, pd))./zq));
g3pred = (2.*3.*f1^2).*(pd*(imag(delfstar2layer(3, d1, phi, pd))./zq));
g5pred = (2.*5.*f1^2).*(pd*(imag(delfstar2layer(5, d1, phi, pd))./zq));
g7pred = (2.*7.*f1^2).*(pd*(imag(delfstar2layer(7, d1, phi, pd))./zq));

set(handles.delf1,'string',num2str(round(f1pred)));
set(handles.delf3,'string',num2str(round(f3pred)));
set(handles.delf5,'string',num2str(round(f5pred)));
set(handles.delf7,'string',num2str(round(f7pred)));

set(handles.delg1,'string',num2str(round(g1pred)));
set(handles.delg3,'string',num2str(round(g3pred)));
set(handles.delg5,'string',num2str(round(g5pred)));
set(handles.delg7,'string',num2str(round(g7pred)));

set(handles.d1,'string',num2str(round(d1,4)));
set(handles.d3,'string',num2str(round(d3,4)));
set(handles.d5,'string',num2str(round(d5,4)));
set(handles.d7,'string',num2str(round(d7,4)));

set(handles.pd1,'string',num2str(round(pdecay1,2)));
set(handles.pd3,'string',num2str(round(pdecay3,2)));
set(handles.pd5,'string',num2str(round(pdecay5,2)));
set(handles.pd7,'string',num2str(round(pdecay7,2)));

set(handles.dev3percent,'string',num2str(round(dev3,2)));
set(handles.dev5percent,'string',num2str(round(dev5,2)));

function lrho = lambdarhof(refn, n, grho, phi)
if refn == n
    lrho = 1/(n*5e6)*(grho^0.5)/cosd(phi/2);
else
    lrho = 1/(n*5e6)*((grho*(n^(phi/90))/(refn^(phi/90)))^0.5)/cosd(phi/2);
end


function d1_Callback(hObject, eventdata, handles)
% hObject    handle to d1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d1 as text
%        str2double(get(hObject,'String')) returns contents of d1 as a double


% --- Executes during object creation, after setting all properties.
function d1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d3_Callback(hObject, eventdata, handles)
% hObject    handle to d3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d3 as text
%        str2double(get(hObject,'String')) returns contents of d3 as a double


% --- Executes during object creation, after setting all properties.
function d3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d5_Callback(hObject, eventdata, handles)
% hObject    handle to d5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d5 as text
%        str2double(get(hObject,'String')) returns contents of d5 as a double


% --- Executes during object creation, after setting all properties.
function d5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pd1_Callback(hObject, eventdata, handles)
% hObject    handle to pd1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pd1 as text
%        str2double(get(hObject,'String')) returns contents of pd1 as a double


% --- Executes during object creation, after setting all properties.
function pd1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pd1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pd3_Callback(hObject, eventdata, handles)
% hObject    handle to pd3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pd3 as text
%        str2double(get(hObject,'String')) returns contents of pd3 as a double


% --- Executes during object creation, after setting all properties.
function pd3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pd3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pd5_Callback(hObject, eventdata, handles)
% hObject    handle to pd5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pd5 as text
%        str2double(get(hObject,'String')) returns contents of pd5 as a double


% --- Executes during object creation, after setting all properties.
function pd5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pd5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d7_Callback(hObject, eventdata, handles)
% hObject    handle to d7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d7 as text
%        str2double(get(hObject,'String')) returns contents of d7 as a double


% --- Executes during object creation, after setting all properties.
function d7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pd7_Callback(hObject, eventdata, handles)
% hObject    handle to pd7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pd7 as text
%        str2double(get(hObject,'String')) returns contents of pd7 as a double


% --- Executes during object creation, after setting all properties.
function pd7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pd7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delf7_Callback(hObject, eventdata, handles)
% hObject    handle to delf7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delf7 as text
%        str2double(get(hObject,'String')) returns contents of delf7 as a double


% --- Executes during object creation, after setting all properties.
function delf7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delf7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delg7_Callback(hObject, eventdata, handles)
% hObject    handle to delg7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delg7 as text
%        str2double(get(hObject,'String')) returns contents of delg7 as a double


% --- Executes during object creation, after setting all properties.
function delg7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delg7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dev3percent_Callback(hObject, eventdata, handles)
% hObject    handle to dev3percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dev3percent as text
%        str2double(get(hObject,'String')) returns contents of dev3percent as a double


% --- Executes during object creation, after setting all properties.
function dev3percent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dev3percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dev5percent_Callback(hObject, eventdata, handles)
% hObject    handle to dev5percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dev5percent as text
%        str2double(get(hObject,'String')) returns contents of dev5percent as a double


% --- Executes during object creation, after setting all properties.
function dev5percent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dev5percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
