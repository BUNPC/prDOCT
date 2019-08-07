%% Phase Resolved Doppler OCT for capillary axial blood flow velocity imaging, GPU-based
% input: 
    % 1D array spectrum, nK*nXrpt*nX*nY, data format: ASCII int16
        % nK: spectrum pixel (camera elements); nXrpt: Ascan repeat;
        % nX: number of Ascans per Bscan; nY: number of Bscans for the whole volum
        % NOTE: the raw data for the whole volume is usually very large, it's recommended to process chunk by chunk
    % PRSinfo: processing information
    % PRSinfo.FWHM: Full width at Half Maxim, Amplitude, [transverse, axial], m
    % PRSinfo.fAline: DAQ Aline rate, Hz
    % PRSinfo.Lam: [light source center, wavelength bandwidth], m
    % PRSinfo.Dim: [nz,nx,nyPchk,nTau]
    % PRSinfo.intDk: OCT lambda to k space interpolation factor (calibration is required)
% subFunctions:
    % function [Dim, fNameBase, fIndex]=GetNameInfoRaw(filename0)
    % function DAT= ReadDat_int16(filePath, Dim, iseg, ARpt_extract,RptBscan) 
    % function RR = DAT2RR_GPU(Dat, intpDk)
    % function [Vz, aRBC]=RR2Vz_GPU(RR, PRSinfo)
% output:
    % Vz, mm/s, [nz,nx,ny]
close all; clear;clc 
%% Default parameter
DefPath  = 'H:\BU\PROJ - Capillary Transit Time\20190716-PMP1\ROI3-DLSOCT-Intralipid';
[filename,datapath]=uigetfile(fullfile(DefPath,'*.dat'),'select file');
dx=1.5; % x pixel size, um
dz=2.7; % z pixel size, um
Crange=[-3 3]; % velocity cAxis range
%% add MATLAB functions' path
addpath('.\SubFunctions') % sub functions
%% 0. data information
[Dim, fNameBase,fIndex]=GetNameInfoRaw(filename);
N_YChks=Dim.ny/10;   % default number of subYseg
%% 1. data process options
clear inputNum
prompt1={['Num Ysegs (nY=',num2str(Dim.ny),', nX=',num2str(Dim.nx),', NxRpt=',num2str(Dim.nxRpt),', NyRpt=',num2str(Dim.nyRpt),')'],...
    'RptA_Start (nARpt Process)','RptA_Inerval (nARpt Process)','RptA_n (nARpt Process)',...
    'N_Aline(DAQ)','N_Bscan(DAQ)','N_xRpt(DAQ)','N_yRpt(DAQ)',...
    'intDk','Aline rate (kHz)'};
inputNum1=inputdlg(prompt1,'', 1,{num2str(N_YChks),...
    '1','1',num2str(max(Dim.nxRpt,Dim.nyRpt)), ...
    num2str(Dim.nx),num2str(Dim.ny),num2str(Dim.nxRpt),num2str(Dim.nyRpt),...
    '-0.43','46'});
N_YChks=str2num(inputNum1{1});    % number of chunks
% extract repeated alines from DAQ nxRpt, for repeated Ascan only
RptA_start=str2num(inputNum1{2});    % selecte start Aline repeat 
RptA_Interval=str2num(inputNum1{3});    % interval 
RptA_n=str2num(inputNum1{4});        % total number of extracted Alines
% DAQ info
Num_Aline=str2num(inputNum1{5});  % number of Aline
Num_Bscan=str2num(inputNum1{6});  % number of Bscans
n_xRpt=str2num(inputNum1{7});     % number of Ascan repeat
n_yRpt=str2num(inputNum1{8});     % number of Bscan repeat
intDk=str2num(inputNum1{9});  % optional
PRSinfo.fAline=str2num(inputNum1{10})*1e3;  % Hz,
PRSinfo.Lam=[1310 170]*1e-9;  % light source [centerWavelength Bandwidth], m
PRSinfo.HPknl=HP_filter_kernel(2.5*1e-3,1/PRSinfo.fAline); % FC=100HZ -> Omega=5*1e-3; FC=200HZ -> Omega=2.5*1e-3; FC=500HZ -> Omega=1*1e-3; 
ARpt_extract=[RptA_start,RptA_Interval,RptA_n];
%% 2. Select focus range %%%%%%
filePath=[datapath,filename];
disp(['Loading data to select the field of focus... ', num2str(floor(N_YChks/2)), ', ',datestr(now,'DD-HH:MM:SS')]);
DimChk=Dim;  DimChk.ny=1;
DAT = ReadDat_int16(filePath, DimChk, floor(N_YChks/2),ARpt_extract); % NII_ori: nk_Nx_ny,Nx=nt*nx;  floor(N_kfile/2)
RR=DAT2RR_GPU(DAT,intDk);  % GPU
clear DAT
fig=figure;
imagesc(abs((squeeze(max(RR(:,:,:),[],3))))); caxis([0 5])
xlabel('Y');ylabel('Z');ylim([1 300]);title('MIP along X')
% caxis([0 5])
disp(['Select brain surface and stack start layer in figure']);
[XY_surf, Z_surf]=ginput(3);
close(fig);
prompt2={'Surface','Start Z_seg','End Z_seg'};
inputZseg=inputdlg(prompt2,'Z Segment parameter', 1,{num2str(floor(Z_surf(1))),num2str(floor(Z_surf(2))),num2str(floor(Z_surf(3)))});
z_seg_surf=str2num(inputZseg{1});
z_seg0=str2num(inputZseg{2});  % number of segments
LengthZ=str2num(inputZseg{3})-str2num(inputZseg{2});
%% 3. Load and process data chunk by chunk
zRange=[z_seg0,z_seg0+LengthZ-1];
ARpt_extract=[RptA_start,RptA_Interval,RptA_n];
nyPerChk=floor(Dim.ny/N_YChks); % number of Bscans per ikfile
DimChk=Dim; DimChk.ny=nyPerChk;
for iChk=1:N_YChks
    disp(['Start loading ith Chunk - ', num2str(iChk), ', ',datestr(now,'DD-HH:MM:SS')]);
    tic;
    [DAT] = ReadDat_int16(filePath, DimChk, iChk, ARpt_extract); %, load the iChk chunk of data, nk_Nx_ny,Nx=nxRpt*nx 
    disp(['DAT2RR...',datestr(now,'DD-HH:MM:SS')]);
    RR0 = DAT2RR_GPU(DAT, intDk);     % GPU-based, process raw spectrum data to spatial reflectivity data 
    RR=permute(reshape(RR0(zRange(1):zRange(2),:,:),[LengthZ,DimChk.nxRpt,DimChk.nx,DimChk.nyRpt,DimChk.ny]),[1 3 5 2 4]); % reshape RR from [nz,Nx,ny] to [nz,nx,ny,nxRpt,nyRpt], Nx=nx*nxRpt
    [nz,nx,nyChk,nxRpt]=size(RR);
    disp(['RR2Vz...',datestr(now,'DD-HH:MM:SS')]);
    Vz(:,:,(iChk-1)*nyPerChk+1:iChk*nyPerChk)=RR2Vz_GPU(RR, PRSinfo)*1e3; toc  % GPU-based prDOCT processing, ~25 times faster than CPU, mm/s
end
save([datapath,'prVz',filename(4:end-4),'.mat'],'-v7.3','Vz')
%%%%%%%%%%%% Vz Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap]=Colormaps_DLSOCT;
Vz3D=imgaussfilt3(Vz,0.5);
zStart=1; zEnd=LengthZ;
figure;
imagesc(squeeze(max(abs(Vz3D(zStart:zEnd,:,:)),[],1)).*sign(squeeze(mean(Vz3D(zStart:zEnd,:,:),1)))); 
colormap(Vzcmap); caxis(Crange); colorbar
axis equal tight
saveas(gcf,[datapath,filename(1:end-4),'.fig']);
saveas(gcf,[datapath,filename(1:end-4),'.png']);