function varargout = GUI_PLOT_prDOCT(varargin)
% GUI_PLOT_prDOCT MATLAB code for GUI_PLOT_prDOCT.fig
%      GUI_PLOT_prDOCT, by itself, creates a new GUI_PLOT_prDOCT or raises the existing
%      singleton*.
%
%      H = GUI_PLOT_prDOCT returns the handle to a new GUI_PLOT_prDOCT or the handle to
%      the existing singleton*.
%
%      GUI_PLOT_prDOCT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PLOT_prDOCT.M with the given input arguments.
%
%      GUI_PLOT_prDOCT('Property','Value',...) creates a new GUI_PLOT_prDOCT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_PLOT_prDOCT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_PLOT_prDOCT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_PLOT_prDOCT

% Last Modified by GUIDE v2.5 17-Jul-2019 14:20:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_PLOT_prDOCT_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_PLOT_prDOCT_OutputFcn, ...
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


% --- Executes just before GUI_PLOT_prDOCT is made visible.
function GUI_PLOT_prDOCT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_PLOT_prDOCT (see VARARGIN)
%% add MATLAB functions' path
handles.CodePath=pwd;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
handles.defpath='H:';

handles.startZ=1;
handles.stackZ=100;
handles.ShowSide=0;
handles.cRange=2;
% Choose default command line output for GUI_PLOT_prDOCT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_PLOT_prDOCT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_PLOT_prDOCT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_LoadData.
function btn_LoadData_Callback(hObject, eventdata, handles)
% hObject    handle to btn_LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
%% select file path %%%%%%%%%
defaultpath=handles.defpath;
[filename,datapath]=uigetfile(defaultpath);
handles.defpath=datapath;
handles.filename=filename;
guidata(hObject, handles);
%% load data %%%%%%%%%%
%%%%% input number of sub GG to be loaded %%%%%%%%
lding=msgbox(['Loading data...  ',datestr(now,'DD:HH:MM')]);
V=LoadMAT(datapath,filename);
lded=msgbox(['Data loaded. ',datestr(now,'DD:HH:MM')]);
pause(1);
delete(lding); delete(lded);
[nz,nx,ny]=size(V);
%%%%%%%%%%%%%
prompt={'dX (um): ', 'dY(um): ', 'dZ size (um)'};
name='Enter Imaging info';
defaultvalue={'1.5','1.5','2.9'};
dXYZinput=inputdlg(prompt,name, 1, defaultvalue);
handles.Xcoor=[1:nx]*str2num(dXYZinput{1});
handles.Ycoor=[1:ny]*str2num(dXYZinput{2});
handles.Zcoor=[1:nz]*str2num(dXYZinput{3});

%% plot velocity enface MIP
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
handles.V=imgaussfilt3(V,0.5);   
imagesc(squeeze(max(abs(handles.V(:,:,:)),[],1)).*sign(squeeze(mean(handles.V(:,:,:),1)))); 
colormap(Vzcmap); caxis([-2 2]); colorbar
axis equal tight
title('V [mm/s]')

guidata(hObject, handles);


% --- Executes on button press in btn_plot.
function btn_plot_Callback(hObject, eventdata, handles)
% hObject    handle to btn_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[nz,nx,ny]=size(handles.V);

prompt={'SideView(N:0;Y:1)','MIP zStart', 'MIP zStack','cRange'};
name='Enter Imaging info';
defaultvalue={num2str(handles.ShowSide),num2str(handles.startZ),num2str(min(nz,nz-handles.startZ)),num2str(handles.cRange)};
dXYZinput=inputdlg(prompt,name, 1, defaultvalue);
handles.ShowSide=str2num(dXYZinput{1});
zStart=str2num(dXYZinput{2});
zStack=str2num(dXYZinput{3});
handles.cRange=str2num(dXYZinput{4});

handles.startZ=zStart;
handles.endZ=min((zStart+zStack)-1,nz);
guidata(hObject, handles);  

%% PLOT
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
if handles.ShowSide==1 % plot SideView figures
    axes(handles.axes1);
    handles.slt(1)=str2num(get(handles.zStart,'string'));
    [handles.slt(3), handles.slt(2)]=ginput(1); % [y x]
    handles.slt=round(handles.slt);
    handles.Show=['3DSP-','[',num2str(handles.slt), ']']; % 3D single plane show
    
    handles.fig=figure;
    set(handles.fig,'Position',[400 500 1400 800])
    subplot(2,2,1) % XY single plan enface view 
    Vxy=(squeeze(max(abs(handles.V(handles.slt(1),:,:)),[],1)).*sign(squeeze(mean(handles.V(handles.slt(1),:,:),1)))); 
    Vxy(handles.slt(2)+[0 1],:)=max(Vxy(:));
    Vxy(:,handles.slt(3)+[0 1])=max(Vxy(:));
    imagesc(handles.Xcoor,handles.Ycoor,Vxy);
    colorbar; colormap(Vzcmap);caxis([-handles.cRange handles.cRange]);
    axis equal; axis tight;
    title(['V-XY, iz=',num2str(handles.slt(1))])
    xlabel('X [um]'); ylabel('Y [um]')
    
    subplot(2,2,2) % xz side view
    Vxz=abs(squeeze(handles.V(:,handles.slt(2),:))).*sign(squeeze(handles.V(:,handles.slt(2),:)));
    Vxz(handles.slt(1),:)=max(Vxz(:));
    imagesc(handles.Ycoor,handles.Zcoor,Vxz);
    colorbar; colormap(Vzcmap);caxis([-handles.cRange handles.cRange]);
    axis equal; axis tight;
    title(['V-XZ, iY=',num2str(handles.slt(2))])
    xlabel('X [um]'); ylabel('Z [um]')
    
    subplot(2,2,3) % XY enface view MIP
    Vsign=sign(squeeze(mean(handles.V(handles.startZ:handles.endZ,:,:),1)));
    Vsign(Vsign~=-1)=1;
    VMIP=(squeeze(max(abs(handles.V(handles.startZ:handles.endZ,:,:)),[],1)).*Vsign); 
    VMIP(handles.slt(2)+[0 1],:)=max(VMIP(:));
    VMIP(:,handles.slt(3)+[0 1])=max(VMIP(:));
    imagesc(handles.Xcoor,handles.Ycoor,VMIP);
    colorbar; colormap(Vzcmap);caxis([-handles.cRange handles.cRange]);
    axis equal; axis tight;
    title(['V MIP [',num2str([handles.startZ handles.endZ]), '] [mm/s]'])
    xlabel('X [um]'); ylabel('Y [um]')
    
    subplot(2,2,4) %yz side view
    Vyz=abs(squeeze(handles.V(:,:,handles.slt(3)))).*sign(squeeze(handles.V(:,:,handles.slt(3))));
    Vyz(handles.slt(1),:)=max(Vyz(:));
    imagesc(handles.Xcoor,handles.Zcoor,Vyz);
    colorbar; colormap(Vzcmap);caxis([-handles.cRange handles.cRange]);
    axis equal; axis tight;
    title(['V-YZ, iX=',num2str(handles.slt(3))])
    xlabel('Y [um]'); ylabel('Z [um]')
else %% plot enface MIP only
    handles.Show=['MIP-[',num2str([handles.startZ handles.endZ]), ']']; % en face MIP show
    handles.fig=figure;
    Vsign=sign(squeeze(mean(handles.V(handles.startZ:handles.endZ,:,:),1)));
    Vsign(Vsign~=-1)=1;
    VMIP=(squeeze(max(abs(handles.V(handles.startZ:handles.endZ,:,:)),[],1)).*Vsign); 
    imagesc(handles.Xcoor,handles.Ycoor,VMIP);
    colorbar; colormap(Vzcmap);caxis([-handles.cRange handles.cRange]);
    axis equal; axis tight;
    title(['V MIP [',num2str([handles.startZ handles.endZ]), '] [mm/s]'])
    xlabel('X [um]'); ylabel('Y [um]')
end
guidata(hObject, handles);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
clc
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
zStart=str2num(get(handles.zStart,'string'));
zStack=str2num(get(handles.zStack,'string'));

[nz,nx,ny] = size(handles.V);
set(hObject,'SliderStep',[1/(nz-1), 3/(nz-1)])
set(hObject,'Max',nz)
zStart=nz-min(round(get(hObject,'Value')),nz-1);
set(handles.zStart,'string',zStart);

VMIP=(squeeze(max(abs(handles.V(zStart:min(zStart+zStack-1,nz),:,:)),[],1)).*sign(squeeze(mean(handles.V(zStart:min(zStart+zStack-1,nz),:,:),1))));
axes(handles.axes1)
imagesc(VMIP);
colorbar; colormap(Vzcmap);caxis([-2 2]);
axis equal; axis tight;
title(['V-XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
xlabel('X [pix]'); ylabel('Y [pix]')
guidata(hObject, handles);

function zStart_Callback(hObject, eventdata, handles)
% hObject    handle to zStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zStart as text
%        str2double(get(hObject,'String')) returns contents of zStart as a double
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
zStart=str2num(get(handles.zStart,'string'));
zStack=str2num(get(handles.zStack,'string'));
[nz,nx,ny]=size(handles.V);

VMIP=(squeeze(max(abs(handles.V(zStart:min(zStart+zStack-1,nz),:,:)),[],1)).*sign(squeeze(mean(handles.V(zStart:min(zStart+zStack-1,nz),:,:),1))));
axes(handles.axes1)
imagesc(VMIP);
colorbar; colormap(Vzcmap);caxis([-2 2]);
axis equal; axis tight;
title(['V-XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
xlabel('X [pix]'); ylabel('Y [pix]')

function zStack_Callback(hObject, eventdata, handles)
% hObject    handle to zStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zStack as text
%        str2double(get(hObject,'String')) returns contents of zStack as a double
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])

[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
zStart=str2num(get(handles.zStart,'string'));
zStack=str2num(get(handles.zStack,'string'));
[nz,nx,ny]=size(handles.V);

VMIP=(squeeze(max(abs(handles.V(zStart:min(zStart+zStack-1,nz),:,:)),[],1)).*sign(squeeze(mean(handles.V(zStart:min(zStart+zStack-1,nz),:,:),1))));
axes(handles.axes1)
imagesc(VMIP);
colorbar; colormap(Vzcmap);caxis([-2 2]);
axis equal; axis tight;
title(['V-XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
xlabel('X [pix]'); ylabel('Y [pix]')

function btn_Save_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
%% select file path %%%%%%%%%
defaultpath=handles.defpath;
disp('Saving data......')
saveas(handles.fig,[defaultpath,handles.DataSlt,'-', handles.Show, '.fig']);
saveas(handles.fig,[defaultpath,handles.DataSlt,'-', handles.Show, '.jpg']);
disp('Data saved!')



% --- Executes on button press in btn_reset.
function btn_reset_Callback(hObject, eventdata, handles)
% hObject    handle to btn_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
handles.defpath='F:';

handles.startZ=1;
handles.stackZ=100;
% Choose default command line output for GUI_PLOT_prDOCT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function zStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes during object creation, after setting all properties.
function zStack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
