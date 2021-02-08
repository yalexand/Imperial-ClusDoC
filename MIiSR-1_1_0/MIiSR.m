function varargout = MIiSR(varargin)
% MIISR MATLAB code for MIiSR.fig
%      MIISR, by itself, creates a new MIISR or raises the existing
%      singleton*.
%
%      H = MIISR returns the handle to a new MIISR or the handle to
%      the existing singleton*.
%
%      MIISR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MIISR.M with the given input arguments.
%
%      MIISR('Property','Value',...) creates a new MIISR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MIiSR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MIiSR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
%==========================================================================
%                             CITATION
%
% This script is provided as a supplemental material in:
%
%   Fabiana A. Caetano, Brennan S. Dirk, Joshua H.K. Tam, P. Craig 
%       Cavanagh, Maria Goiko, Stephen S.G. Ferguson, Stephen H. Pasternak,
%       Jimmy D. Dikeakos, John R. de Bruyn, Bryan Heit. MIiSR: Analysis of 
%		Molecular Interactions in Super-Resolution Imaging Enables the Study 
%		of Protein Interactions, Dynamics and Formation of Multi-protein 
%		Structures. 2015. PLoS Computational Biology
% 
% Please reference this paper in any publications which use this script for 
% analysis.
%==========================================================================

% Edit the above text to modify the response to help MIiSR

% Last Modified by GUIDE v2.5 01-Oct-2015 14:00:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MIiSR_OpeningFcn, ...
                   'gui_OutputFcn',  @MIiSR_OutputFcn, ...
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


% --- Executes just before MIiSR is made visible.
function MIiSR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MIiSR (see VARARGIN)

% Choose default command line output for MIiSR
handles.output = hObject;
handles.Nqueue = 0; %number of entries in the queue
handles.pwd = pwd; %home path
axes(handles.imDisplay);
imshow([]);

%check for active parallel pool
b = gcp('nocreate');

if ~isempty(b)
    nCores = b.NumWorkers;
    set (handles.nCoreText, 'String', num2str(nCores));
    set (handles.stopPool, 'Enable', 'on');
    set (handles.startPool, 'Enable', 'off');
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MIiSR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MIiSR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in addQueue.
function addQueue_Callback(hObject, eventdata, handles)
% hObject    handle to addQueue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.addCh1, 'Enable', 'off');
set(handles.addCh2, 'Enable', 'off');
set(handles.addCh3, 'Enable', 'off');
set(handles.Ch1Name, 'Enable', 'off');
set(handles.Ch2Name, 'Enable', 'off');
set(handles.Ch3Name, 'Enable', 'off');
set(handles.loadIm, 'Enable', 'off');
set(handles.newROI, 'Enable', 'off');
set(handles.addQueue, 'Enable', 'off');
set(handles.startButton, 'Enable', 'off');
set(handles.clearImage, 'Enable', 'off');

handles.Nqueue = handles.Nqueue + 1; %set position in queue

% --------generate structure, add to handle.queue
%files/filenames
tPath = get(handles.Ch1Path, 'String');
s =  strfind(tPath, filesep);
s = s(end);
tPath(s:end) = [];
handles.queue(handles.Nqueue).path = tPath;

handles.queue(handles.Nqueue).nCh = handles.nCh;
handles.queue(handles.Nqueue).Ch1Name = get(handles.Ch1Name, 'String');
handles.queue(handles.Nqueue).Ch1Path = get(handles.Ch1Path, 'String');

if handles.nCh == 2
    handles.queue(handles.Nqueue).Ch2Name = get(handles.Ch2Name, 'String');
    handles.queue(handles.Nqueue).Ch2Path = get(handles.Ch2Path, 'String');
    handles.queue(handles.Nqueue).Ch3Name = [];
    handles.queue(handles.Nqueue).Ch3Path = [];
elseif handles.nCh == 3
    handles.queue(handles.Nqueue).Ch2Name = get(handles.Ch2Name, 'String');
    handles.queue(handles.Nqueue).Ch2Path = get(handles.Ch2Path, 'String');
    handles.queue(handles.Nqueue).Ch3Name = get(handles.Ch3Name, 'String');
    handles.queue(handles.Nqueue).Ch3Path = get(handles.Ch3Path, 'String');
else
    handles.queue(handles.Nqueue).Ch2Name = [];
    handles.queue(handles.Nqueue).Ch2Path = [];
    handles.queue(handles.Nqueue).Ch3Name = [];
    handles.queue(handles.Nqueue).Ch3Path = [];
end

%store values for filtering
handles.queue(handles.Nqueue).denistyFilter = get(handles.denistyFilter, 'Value');
handles.queue(handles.Nqueue).densitySD1 = str2double(get(handles.densitySD1, 'String'));
handles.queue(handles.Nqueue).densitySD2 = str2double(get(handles.densitySD2, 'String'));
handles.queue(handles.Nqueue).densitySD3 = str2double(get(handles.densitySD3, 'String'));
handles.queue(handles.Nqueue).densityDist = str2double(get(handles.densityDist, 'String'));

%store values for SAA
handles.queue(handles.Nqueue).SAAcheck = get(handles.SAAcheck, 'Value');
handles.queue(handles.Nqueue).SAAdist = get(handles.SAAdist, 'String');
handles.queue(handles.Nqueue).SAArand = get(handles.SAArand, 'String');
handles.queue(handles.Nqueue).SAAaF = get(handles.SAAaF, 'String');
handles.queue(handles.Nqueue).polyRand = get(handles.polyRand, 'Value');
handles.queue(handles.Nqueue).calcCDC = get(handles.calcCDC, 'Value');
handles.queue(handles.Nqueue).CDC_SD = get(handles.CDC_SD, 'String');
handles.queue(handles.Nqueue).SAAcutoff = get(handles.SAAcutoff, 'String');
handles.queue(handles.Nqueue).SAAuserCDC = get(handles.SAAuserCDC, 'Value');
handles.queue(handles.Nqueue).SAAuserCDCvalue = get(handles.SAAuserCDCvalue, 'String');

%store values for spatial statistics
handles.queue(handles.Nqueue).RDFcheck = get(handles.RDFcheck, 'Value');
handles.queue(handles.Nqueue).spatialDist = get(handles.spatialDist, 'String');
handles.queue(handles.Nqueue).spatialCh1 = get(handles.spatialCh1, 'Value');
handles.queue(handles.Nqueue).spatialCh2 = get(handles.spatialCh2, 'Value');

%store values for clustering analysis
handles.queue(handles.Nqueue).DBSCANcheck = get(handles.DBSCANcheck, 'Value');
handles.queue(handles.Nqueue).DBSCANk = get(handles.DBSCANk, 'String');
handles.queue(handles.Nqueue).DBSCANe = get(handles.DBSCANe, 'String');
handles.queue(handles.Nqueue).DBSCANch = get(handles.DBSCANch, 'Value');
handles.queue(handles.Nqueue).OPTICScheck = get(handles.OPTICScheck, 'Value');
handles.queue(handles.Nqueue).OPTICSk = get(handles.OPTICSk, 'String');
handles.queue(handles.Nqueue).OPTICSch = get(handles.OPTICSch, 'Value');
handles.queue(handles.Nqueue).OPTICShier = get(handles.OPTICShier, 'Value');
handles.queue(handles.Nqueue).hierRatio = get(handles.hierRatio, 'String');

%determine file name, write info to queue list
handles.queue(handles.Nqueue).cropCoords = handles.ROIposition;
qList = get(handles.queueList, 'String'); %get current display queue
if handles.nCh == 1
    fName = handles.queue(handles.Nqueue).Ch1Name;
elseif handles.nCh == 2
    fName = [handles.queue(handles.Nqueue).Ch1Name '_' handles.queue(handles.Nqueue).Ch2Name];
else
    fName = [handles.queue(handles.Nqueue).Ch1Name '_' handles.queue(handles.Nqueue).Ch2Name '_' handles.queue(handles.Nqueue).Ch3Name];
end
tName = cellstr([fName ' - [' num2str(handles.ROIposition(1)) ', ' num2str(handles.ROIposition(2)) ', ' num2str(handles.ROIposition(3)) ', ' num2str(handles.ROIposition(4)) ']']);
%check to see if filename exists
if size(qList,1) >= 1 && ~isempty(qList)
    nameExists = strfind(qList, char(tName));
    nameExists = cell2mat(nameExists);
    % add suffix if files with same name exists
    if ~isempty(nameExists)
        fExt = size(nameExists,1)+1;
        tName = char(tName);
        tName = [tName ' - ' num2str(fExt)];
    end
end

handles.queue(handles.Nqueue).fName = tName;
qList = [qList; tName];
set(handles.queueList, 'String', qList);

%reset for next ROI
set(handles.addCh1, 'Enable', 'on');
set(handles.Ch1Name, 'Enable', 'on');

if handles.nCh > 1
    set(handles.addCh2, 'Enable', 'on');
    set(handles.Ch2Name, 'Enable', 'on');
end
if handles.nCh == 3
    set(handles.addCh3, 'Enable', 'on');
    set(handles.Ch3Name, 'Enable', 'on');
end

set(handles.loadIm, 'Enable', 'on');
set(handles.newROI, 'Enable', 'on');
set(handles.startButton, 'Enable', 'on');
set(handles.addQueue, 'Enable', 'on');
set(handles.clearImage, 'Enable', 'on');
guidata(hObject, handles);

% --- Executes on button press in delQueue.
function delQueue_Callback(hObject, eventdata, handles)
% hObject    handle to delQueue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selLine = get(handles.queueList, 'Value');
qList = get(handles.queueList, 'String'); %get current display queue
qList(selLine) = [];
set(handles.queueList, 'String', qList); %rewrite queue
pause (0.5);
handles.queue(selLine) = []; %remove queue entry
handles.Nqueue = handles.Nqueue - 1; %set position in queue
guidata(hObject, handles);


% --- Executes on selection change in queueList.
function queueList_Callback(hObject, eventdata, handles)
% hObject    handle to queueList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns queueList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from queueList


% --- Executes during object creation, after setting all properties.
function queueList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to queueList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
HelpFile = which('MIiSR.m');
s=strfind(HelpFile,'MIiSR.m');
HelpFile(s:end)=[];
HelpFile = [HelpFile 'help' filesep 'MIiSRHelp.html'];
web (HelpFile);

% --- Executes on button press in addCh1.
function addCh1_Callback(hObject, eventdata, handles)
% hObject    handle to addCh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[tmpFile, tmpPath] = uigetfile('*.mat', 'Select a .mat file for channel 1'); % get file

if length(tmpFile)>1
    tmpFile = [tmpPath, tmpFile];
    set(handles.Ch1Path, 'String', tmpFile);
    set(handles.addCh2, 'Enable', 'on');
    set(handles.Ch2Name, 'Enable', 'on');
    cd (tmpPath);
    guidata(hObject, handles);
else
    %clear & inactive buttons
    set(handles.Ch1Path, 'String', []);
    set(handles.Ch2Path, 'String', []);
    set(handles.Ch2Name, 'Enable', 'off');
    set(handles.addCh2, 'Enable', 'off');
    set(handles.Ch3Path, 'String', []);
    set(handles.Ch3Name, 'Enable', 'off');
    set(handles.addCh3, 'Enable', 'off');
    guidata(hObject, handles);
end

% --- Executes on button press in addCh2.
function addCh2_Callback(hObject, eventdata, handles)
% hObject    handle to addCh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[tmpFile, tmpPath] = uigetfile('*.mat', 'Select a .mat file for channel 2'); % get file

if length(tmpFile)>1
    tmpFile = [tmpPath, tmpFile];
    set(handles.Ch2Path, 'String', tmpFile);
    set(handles.addCh3, 'Enable', 'on');
    set(handles.Ch3Name, 'Enable', 'on');
    cd (tmpPath);
    guidata(hObject, handles);
else
    %clear & inactive buttons
    set(handles.Ch2Path, 'String', []);
    set(handles.Ch3Path, 'String', []);
    set(handles.Ch3Name, 'Enable', 'off');
    set(handles.addCh3, 'Enable', 'off');
    guidata(hObject, handles);
end


% --- Executes on button press in addCh3.
function addCh3_Callback(hObject, eventdata, handles)
% hObject    handle to addCh3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[tmpFile, tmpPath] = uigetfile('*.mat', 'Select a .mat file for channel 3'); % get file

if length(tmpFile)>1
    tmpFile = [tmpPath, tmpFile];
    cd (tmpPath);
    set(handles.Ch3Path, 'String', tmpFile);
    guidata(hObject, handles);
else
    %clear & inactive buttons
    set(handles.Ch3Path, 'String', []);
    guidata(hObject, handles);
end


% --- Executes on button press in loadIm.
function loadIm_Callback(hObject, eventdata, handles)
% hObject    handle to loadIm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% load image into temporary file

%inactive file control portion of GUI
set(handles.addCh1, 'Enable', 'off');
set(handles.addCh2, 'Enable', 'off');
set(handles.addCh3, 'Enable', 'off');
set(handles.Ch1Name, 'Enable', 'off');
set(handles.Ch2Name, 'Enable', 'off');
set(handles.Ch3Name, 'Enable', 'off');
set(handles.loadIm, 'Enable', 'off');
set(handles.newROI, 'Enable', 'off');
set(handles.delROI, 'Enable', 'off');
set(handles.clearImage, 'Enable', 'off');
set(handles.startButton, 'Enable', 'off');


%load position files
Im1 = get(handles.Ch1Path, 'String');
if ~isempty(Im1) %file selected
    %move to directory of channel 1 (all files saved here)
    aDir = char(get(handles.Ch1Path, 'String'));
    s = strfind(aDir, filesep);
    s = s(end);
    aDir(s:end) = [];
    cd (aDir);
    pause (0.5);
    Im1 = load(get(handles.Ch1Path, 'String'));
    handles.nCh = 1;
else
    handles.nCh = 0;
end

Im2 = get(handles.Ch2Path, 'String');
if isempty(Im2) %single-channel image
    Im2 = [];
    Im2.outMat = [0,0];
    Im3 = [];
    Im3.outMat = [0,0];
else
    Im2 = load(get(handles.Ch2Path, 'String')); %2-channel image
    handles.nCh = 2;
    Im3 = get(handles.Ch3Path, 'String');
    if isempty(Im3) 
        Im3 = [];
        Im3.outMat = [0,0];
    else %3-channel image
        Im3 = load(get(handles.Ch3Path, 'String'));
        handles.nCh = 3;
    end
end


if handles.nCh == 1
    %scale
    maxX = ceil(max(Im1.outMat(:,1)));
    maxY = ceil(max(Im1.outMat(:,2)));
    maxX = ceil(maxX/20); %scale to 20nm pixels
    maxY = ceil(maxY/20); %scale to 20nm pixels
    tmpIm = zeros(maxX, maxY); %empty grayscale image
    for jj = 1:length(Im1.outMat)
        tPos(1) = ceil(Im1.outMat(jj,1)/20); %scale to 20nm pixels
        tPos(2) = ceil(Im1.outMat(jj,2)/20); %scale to 20nm pixels
        tmpIm(tPos(1),tPos(2),1) = tmpIm(tPos(1),tPos(2),1) + 1;
    end
    %scale to 255
    maxI = max(tmpIm(:));
    tmpIm = (tmpIm./maxI).*255;
    tmpIm = uint8(tmpIm);
elseif handles.nCh > 1
    %scale
    maxX = max([ceil(max(Im1.outMat(:,1))), ceil(max(Im2.outMat(:,1))), ceil(max(Im3.outMat(:,1)))]);
    maxY = max([ceil(max(Im1.outMat(:,2))), ceil(max(Im2.outMat(:,2))), ceil(max(Im3.outMat(:,2)))]);
    maxX = ceil(maxX/20); %scale to 20nm pixels
    maxY = ceil(maxY/20); %scale to 20nm pixels
    tmpIm = zeros(maxX, maxY,3); %empty RGB image
    for jj = 1:length(Im1.outMat)
        tPos(1) = ceil(Im1.outMat(jj,1)/20); %scale to 20nm pixels; ensure no zero indices
        tPos(2) = ceil(Im1.outMat(jj,2)/20); %scale to 20nm pixels; ensure no zero indices
        tmpIm(tPos(1),tPos(2),1) = tmpIm(tPos(1),tPos(2),1) + 1;
    end
    %scale to 255
    maxI = tmpIm(:,:,1);
    maxI = max(maxI(:));
    tmpIm = (tmpIm./maxI).*255;
    
    for jj = 1:length(Im2.outMat)
        tPos(1) = ceil(Im2.outMat(jj,1)/20); %scale to 20nm pixels
        tPos(2) = ceil(Im2.outMat(jj,2)/20); %scale to 20nm pixels
        tmpIm(tPos(1),tPos(2),2) = tmpIm(tPos(1),tPos(2),2) + 1;
    end
    %scale to 255
    maxI = tmpIm(:,:,1);
    maxI = max(maxI(:));
    tmpIm(:,:,1) = (tmpIm(:,:,1)./maxI).*255;
    maxI = tmpIm(:,:,2);
    maxI = max(maxI(:));
    tmpIm(:,:,2) = (tmpIm(:,:,2)./maxI).*255;
    if handles.nCh == 3
        for jj = 1:length(Im3.outMat)
            tPos(1) = ceil(Im3.outMat(jj,1)/20); %scale to 20nm pixels
            tPos(2) = ceil(Im3.outMat(jj,2)/20); %scale to 20nm pixels
            tmpIm(tPos(1),tPos(2),3) = tmpIm(tPos(1),tPos(2),3) + 1;
        end
        maxI = tmpIm(:,:,3);
        maxI = max(maxI(:));
        tmpIm(:,:,3) = (tmpIm(:,:,3)./maxI).*255;
    end
    tmpIm = uint8(tmpIm); %scale to 8-bit RGB image    
end

if handles.nCh > 1
    imMax = double(tmpIm(:,:,1));
    imMax(imMax==0) = [];
    imMax = quantile(imMax(:),4);
    imMax = imMax(3)/100;
    if imMax > 1
        imMax = 1;
    end
    handles.Image = imadjust(tmpIm, [0 0 0; imMax imMax imMax, []]);
    axes(handles.imDisplay);
    imshow(handles.Image);
elseif handles.nCh == 1
    imMax = double(tmpIm(:,:,1));
    imMax(imMax==0) = [];
    imMax = quantile(imMax(:),4);
    handles.Image = tmpIm;
    axes(handles.imDisplay);
    imshow(handles.Image, [0 imMax(3)]);
end    

%reactive file control portion of GUI, set for new channel
set(handles.addCh1, 'Enable', 'on');
set(handles.Ch1Name, 'Enable', 'on');
if handles.nCh == 2
    set(handles.Ch2Name, 'Enable', 'on');
    set(handles.addCh2, 'Enable', 'on');
elseif handles.nCh ==3
    set(handles.Ch2Name, 'Enable', 'on');
    set(handles.addCh2, 'Enable', 'on');
    set(handles.Ch3Name, 'Enable', 'on');
    set(handles.addCh3, 'Enable', 'on');
end
set(handles.loadIm, 'Enable', 'on');
set(handles.newROI, 'Enable', 'on');
set(handles.startButton, 'Enable', 'on');
set(handles.delROI, 'Enable', 'on');
set(handles.clearImage, 'Enable', 'on');
guidata(hObject, handles);

function Ch1Path_Callback(hObject, eventdata, handles)
% hObject    handle to Ch1Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ch1Path as text
%        str2double(get(hObject,'String')) returns contents of Ch1Path as a double


% --- Executes during object creation, after setting all properties.
function Ch1Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ch1Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ch2Path_Callback(hObject, eventdata, handles)
% hObject    handle to Ch2Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ch2Path as text
%        str2double(get(hObject,'String')) returns contents of Ch2Path as a double


% --- Executes during object creation, after setting all properties.
function Ch2Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ch2Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ch3Path_Callback(hObject, eventdata, handles)
% hObject    handle to Ch3Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ch3Path as text
%        str2double(get(hObject,'String')) returns contents of Ch3Path as a double


% --- Executes during object creation, after setting all properties.
function Ch3Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ch3Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ch1Name_Callback(hObject, eventdata, handles)
% hObject    handle to Ch1Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ch1Name as text
%        str2double(get(hObject,'String')) returns contents of Ch1Name as a double


% --- Executes during object creation, after setting all properties.
function Ch1Name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ch1Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ch2Name_Callback(hObject, eventdata, handles)
% hObject    handle to Ch2Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ch2Name as text
%        str2double(get(hObject,'String')) returns contents of Ch2Name as a double


% --- Executes during object creation, after setting all properties.
function Ch2Name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ch2Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ch3Name_Callback(hObject, eventdata, handles)
% hObject    handle to Ch3Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ch3Name as text
%        str2double(get(hObject,'String')) returns contents of Ch3Name as a double


% --- Executes during object creation, after setting all properties.
function Ch3Name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ch3Name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in localCheck.
function localCheck_Callback(hObject, eventdata, handles)
% hObject    handle to localCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of localCheck
lCheck = get(handles.localCheck, 'Value');
if lCheck
    set (handles.NcoreCheck, 'Value', 0);
    set (handles.nCores, 'String', 'n');
    set (handles.nCores, 'Enable', 'off');
else
    set (handles.NcoreCheck, 'Value', 1);
    set (handles.nCores, 'Enable', 'on');
end
guidata(hObject, handles);

% --- Executes on button press in NcoreCheck.
function NcoreCheck_Callback(hObject, eventdata, handles)
% hObject    handle to NcoreCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NcoreCheck
lCheck = get(handles.NcoreCheck, 'Value');
if lCheck
    set (handles.localCheck, 'Value', 0);
    set (handles.nCores, 'Enable', 'on');
else
    set (handles.localCheck, 'Value', 1);
    set (handles.nCores, 'Enable', 'off');
    set (handles.nCores, 'String', 'n');
end
guidata(hObject, handles);

function nCores_Callback(hObject, eventdata, handles)
% hObject    handle to nCores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nCores as text
%        str2double(get(hObject,'String')) returns contents of nCores as a double


% --- Executes during object creation, after setting all properties.
function nCores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nCores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in startPool.
function startPool_Callback(hObject, eventdata, handles)
% hObject    handle to startPool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% turn off UI
set(findobj(handles.uipanel1, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel2, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel4, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel5, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel6, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel9, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel11, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel7, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel13, '-property', 'Enable'), 'Enable', 'off');
pause (0.5); %pause to allow GUI to update


b = gcp('nocreate');
if ~isempty(b)
    delete (gcp);
end

useLocal = get(handles.localCheck, 'Value');
if useLocal
    set (handles.stopPool, 'Enable', 'off');
    set (handles.startPool, 'Enable', 'off');
    parpool();
    b = gcp('nocreate');
    nCores = b.NumWorkers;
    set (handles.nCoreText, 'String', num2str(nCores));
    set (handles.stopPool, 'Enable', 'on');

else
    set (handles.stopPool, 'Enable', 'off');
    set (handles.startPool, 'Enable', 'off');
    set (handles.nCores, 'Enable', 'off');
    nCores = get(handles.nCores, 'String');
    nCores = ceil(str2num(nCores));
    if isempty(nCores) || nCores < 1
        parpool(); % use local
        b = gcp('nocreate');
        nCores = b.NumWorkers;
    else
        parpool(nCores);
    end
    set (handles.nCoreText, 'String', num2str(nCores));
    set (handles.stopPool, 'Enable', 'on');
    set (handles.nCores, 'Enable', 'on');
end

% turn on UI
set(findobj(handles.uipanel1, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel2, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel4, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel5, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel6, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel9, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel11, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel7, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel13, '-property', 'Enable'), 'Enable', 'on');
set(handles.startPool, 'Enable', 'off');
denistyFilter_Callback(hObject, [], handles);
SAAcheck_Callback(hObject, [], handles);
calcCDC_Callback(hObject, [], handles);
polyRand_Callback(hObject, [], handles);
RDFcheck_Callback(hObject, [], handles);
DBSCANcheck_Callback(hObject, [], handles);
OPTICScheck_Callback(hObject, [], handles);
OPTICShier_Callback(hObject, [], handles);
denistyFilter_Callback(hObject, [], handles);
pause (0.5); %pause to allow GUI to update


% --- Executes on button press in stopPool.
function stopPool_Callback(hObject, eventdata, handles)
% hObject    handle to stopPool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% turn off UI
set(findobj(handles.uipanel1, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel2, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel4, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel5, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel6, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel9, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel11, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel7, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel13, '-property', 'Enable'), 'Enable', 'off');
pause (0.5); %pause to allow GUI to update


b = gcp('nocreate');
if ~isempty(b)
    delete (gcp);
end
set (handles.stopPool, 'Enable', 'off');
set (handles.startPool, 'Enable', 'on');
set (handles.nCoreText, 'String', '[]');
% turn on UI
set(findobj(handles.uipanel1, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel2, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel4, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel5, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel6, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel9, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel11, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel7, '-property', 'Enable'), 'Enable', 'on');
set(findobj(handles.uipanel13, '-property', 'Enable'), 'Enable', 'on');
set(handles.stopPool, 'Enable', 'off');
set(handles.nCores, 'Enable', 'off');
set(handles.nCores, 'String', 'n');
denistyFilter_Callback(hObject, [], handles);
SAAcheck_Callback(hObject, [], handles);
calcCDC_Callback(hObject, [], handles);
polyRand_Callback(hObject, [], handles);
RDFcheck_Callback(hObject, [], handles);
DBSCANcheck_Callback(hObject, [], handles);
OPTICScheck_Callback(hObject, [], handles);
OPTICShier_Callback(hObject, [], handles);

% --- Executes on button press in newROI.
function newROI_Callback(hObject, eventdata, handles)
% hObject    handle to newROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield (handles, 'h')
    delete (handles.h);
end
set (handles.newROI, 'Enable', 'off');
set (handles.addQueue, 'Enable', 'off');
axes(handles.imDisplay);
handles.h = imrect(handles.imDisplay);
    % more modern version
    %handles.h = drawrectangle(handles.imDisplay);
tpos = wait(handles.h);
tpos(tpos<0) = 0; %remove any negative indices
ROIpos(1) = tpos(2);
ROIpos(2) = tpos(2)+tpos(4); %convert xmin/max, ymin/max
ROIpos(3) = tpos(1);
ROIpos(4) = tpos(1)+tpos(3);
handles.ROIposition = ROIpos;
set (handles.newROI, 'Enable', 'on');
set (handles.addQueue, 'Enable', 'on');
guidata(hObject, handles);

pause (0.5); %pause to allow GUI to update

% --- Executes on button press in SAAcheck.
function SAAcheck_Callback(hObject, eventdata, handles)
% hObject    handle to SAAcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SAAcheck
tVal = get(handles.SAAcheck, 'Value');
if tVal
    set (handles.SAAuserCDC, 'Value', 0);
    set (handles.calcCDC, 'Value', 1);
    set (handles.SAAuserCDCvalue, 'Enable', 'off');
    set (handles.CDC_SD, 'Enable', 'on');
    set(handles.SAAcutoff, 'Enable', 'on');
    set (handles.calcCDC, 'Enable', 'on');
    set (handles.SAAuserCDC, 'Enable', 'on');
    set (handles.SAAaF, 'Enable', 'on');
    set (handles.SAAdist, 'Enable', 'on');
    set (handles.SAArand, 'Enable', 'on');
else
    set (handles.SAAuserCDC, 'Value', 0);
    set (handles.SAAuserCDC, 'Enable', 'off');
    set (handles.calcCDC, 'Value', 0);
    set (handles.calcCDC, 'Enable', 'off');
    set (handles.SAAuserCDCvalue, 'Enable', 'off');
    set (handles.CDC_SD, 'Enable', 'off');
    set(handles.SAAcutoff, 'Enable', 'off');
    set (handles.SAAaF, 'Enable', 'off');
    set (handles.SAAdist, 'Enable', 'off');
    set (handles.SAArand, 'Enable', 'off');

end

% --- Executes on button press in calcCDC.
function calcCDC_Callback(hObject, eventdata, handles)
% hObject    handle to calcCDC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcCDC
tVal = get(handles.calcCDC, 'Value');
if tVal
    set (handles.SAAuserCDC, 'Value', 0);
    set (handles.calcCDC, 'Value', 1);
    set (handles.SAAuserCDCvalue, 'Enable', 'off');
    set (handles.CDC_SD, 'Enable', 'on');
    set(handles.SAAcutoff, 'Enable', 'on');
else
    set (handles.SAAuserCDC, 'Value', 1);
    set (handles.calcCDC, 'Value', 0);
    set (handles.SAAuserCDCvalue, 'Enable', 'on');
    set (handles.CDC_SD, 'Enable', 'off');
    set(handles.SAAcutoff, 'Enable', 'off');
end


function CDC_SD_Callback(hObject, eventdata, handles)
% hObject    handle to CDC_SD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CDC_SD as text
%        str2double(get(hObject,'String')) returns contents of CDC_SD as a double


% --- Executes during object creation, after setting all properties.
function CDC_SD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CDC_SD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SAAcutoff_Callback(hObject, eventdata, handles)
% hObject    handle to SAAcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SAAcutoff as text
%        str2double(get(hObject,'String')) returns contents of SAAcutoff as a double


% --- Executes during object creation, after setting all properties.
function SAAcutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAAcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SAAdist_Callback(hObject, eventdata, handles)
% hObject    handle to SAAdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SAAdist as text
%        str2double(get(hObject,'String')) returns contents of SAAdist as a double


% --- Executes during object creation, after setting all properties.
function SAAdist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAAdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SAArand_Callback(hObject, eventdata, handles)
% hObject    handle to SAArand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SAArand as text
%        str2double(get(hObject,'String')) returns contents of SAArand as a double


% --- Executes during object creation, after setting all properties.
function SAArand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAArand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SAAuserCDC.
function SAAuserCDC_Callback(hObject, eventdata, handles)
% hObject    handle to SAAuserCDC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SAAuserCDC
tVal = get(handles.SAAuserCDC, 'Value');
if tVal
    set (handles.SAAuserCDC, 'Value', 1);
    set (handles.calcCDC, 'Value', 0);
    set (handles.SAAuserCDCvalue, 'Enable', 'on');
    set (handles.CDC_SD, 'Enable', 'off');
    set(handles.SAAcutoff, 'Enable', 'off');
else
    set (handles.SAAuserCDC, 'Value', 0);
    set (handles.calcCDC, 'Value', 1);
    set (handles.SAAuserCDCvalue, 'Enable', 'off');
    set (handles.CDC_SD, 'Enable', 'on');
    set(handles.SAAcutoff, 'Enable', 'on');
end

function SAAuserCDCvalue_Callback(hObject, eventdata, handles)
% hObject    handle to SAAuserCDCvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SAAuserCDCvalue as text
%        str2double(get(hObject,'String')) returns contents of SAAuserCDCvalue as a double


% --- Executes during object creation, after setting all properties.
function SAAuserCDCvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAAuserCDCvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RDFcheck.
function RDFcheck_Callback(hObject, eventdata, handles)
% hObject    handle to RDFcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RDFcheck



function spatialDist_Callback(hObject, eventdata, handles)
% hObject    handle to spatialDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spatialDist as text
%        str2double(get(hObject,'String')) returns contents of spatialDist as a double


% --- Executes during object creation, after setting all properties.
function spatialDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spatialDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DBSCANcheck.
function DBSCANcheck_Callback(hObject, eventdata, handles)
% hObject    handle to DBSCANcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DBSCANcheck
tVal = get(handles.DBSCANcheck, 'Value');
if tVal
    set(handles.DBSCANk, 'Enable', 'on');
    set(handles.DBSCANe, 'Enable', 'on');
    set(handles.DBSCANch, 'Enable', 'on');
else
    set(handles.DBSCANk, 'Enable', 'off');
    set(handles.DBSCANe, 'Enable', 'off');
    set(handles.DBSCANch, 'Enable', 'off');
end    

% --- Executes on button press in OPTICScheck.
function OPTICScheck_Callback(hObject, eventdata, handles)
% hObject    handle to OPTICScheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OPTICScheck
tVal = get(handles.OPTICScheck, 'Value');
if tVal
    set(handles.OPTICSk, 'Enable', 'on');
    set(handles.OPTICSch, 'Enable', 'on');
    set(handles.OPTICShier, 'Enable', 'on');
    OPTICShier_Callback(hObject, [], handles);
else
    set(handles.OPTICSk, 'Enable', 'off');
    set(handles.OPTICSch, 'Enable', 'off');
    set(handles.OPTICShier, 'Enable', 'off');
    set(handles.hierRatio, 'Enable', 'off');
end   


function DBSCANk_Callback(hObject, eventdata, handles)
% hObject    handle to DBSCANk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DBSCANk as text
%        str2double(get(hObject,'String')) returns contents of DBSCANk as a double


% --- Executes during object creation, after setting all properties.
function DBSCANk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DBSCANk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DBSCANe_Callback(hObject, eventdata, handles)
% hObject    handle to DBSCANe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DBSCANe as text
%        str2double(get(hObject,'String')) returns contents of DBSCANe as a double


% --- Executes during object creation, after setting all properties.
function DBSCANe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DBSCANe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in DBSCANch.
function DBSCANch_Callback(hObject, eventdata, handles)
% hObject    handle to DBSCANch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DBSCANch contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DBSCANch


% --- Executes during object creation, after setting all properties.
function DBSCANch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DBSCANch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function OPTICSk_Callback(hObject, eventdata, handles)
% hObject    handle to OPTICSk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OPTICSk as text
%        str2double(get(hObject,'String')) returns contents of OPTICSk as a double


% --- Executes during object creation, after setting all properties.
function OPTICSk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OPTICSk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in OPTICSch.
function OPTICSch_Callback(hObject, eventdata, handles)
% hObject    handle to OPTICSch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns OPTICSch contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OPTICSch


% --- Executes during object creation, after setting all properties.
function OPTICSch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OPTICSch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OPTICShier.
function OPTICShier_Callback(hObject, eventdata, handles)
% hObject    handle to OPTICShier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OPTICShier
tVal = get(handles.OPTICShier, 'Value');
if tVal
    set(handles.hierRatio, 'Enable', 'on');
else
    set(handles.hierRatio, 'Enable', 'off');
end

% --- Executes on button press in startButton.
function startButton_Callback(hObject, eventdata, handles)
% hObject    handle to startButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% turn off UI
set(findobj(handles.uipanel1, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel2, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel4, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel5, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel6, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel9, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel11, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel7, '-property', 'Enable'), 'Enable', 'off');
set(findobj(handles.uipanel13, '-property', 'Enable'), 'Enable', 'off');
pause (0.5); %pause to allow GUI to update

%start parallel processing, if not already done
b = gcp('nocreate');
if isempty(b)
    startPool_Callback(hObject, [], handles); %start pool if not done already
    pause (0.5);
end

% loop through queue
for ii=1:handles.Nqueue
    display (' ');
    display (['Analyzing ROI ' num2str(ii) ' of ' num2str(handles.Nqueue) '...']);
    a = handles.queue(ii); %temp structure
    MIiSRdata.queueInfo = handles.queue;
    
    %update queue list
    qList = get(handles.queueList, 'String'); %get current display queue
    tName = char(qList(1));
    tName = cellstr(['Processing: ' tName]);
    qList(1) = tName;
    set(handles.queueList, 'String', qList);
    pause (0.5); %provide time for GUI to update
    
    % --------------------   Load & Crop Images
    cd (a.path); %move to base directory
    tName = char(a.fName);
    s = strfind(tName, '[');
    s = s(1)-3; %find divider in file name and crop area
    Conditions.fName = tName(1:s);
    mkdir (tName); %make folder for files
    cd (tName); %move to folder, to save all files.
    
    Conditions.ImScale = 20; %nm/pixel in TIFF images
    Conditions.saveTIFF = 1; %save TIFF files
	Conditions.saveUncropped = 0; %save uncropped image
    Conditions.minPhotons = 0; %filter for minimal number of photons, 0 = no filter
    Conditions.minPrecission = 0; %filter for minimal precision (in nm), 0 = no filter
    Conditions.savePosFiles = 1; %save position (.mat) files
    Conditions.saveSingleCh = 0; %save single channels
	Conditions.Ch1Label = char(a.Ch1Name);
    Conditions.Ch2Label = char(a.Ch2Name);
    Conditions.Ch3Label = char(a.Ch3Name);
    
    % generate file processing list
    if a.nCh == 1
        fList = {a.Ch1Path};
    elseif a.nCh == 2

        fList = {a.Ch1Path; a.Ch2Path};
    else
        fList = {a.Ch1Path; a.Ch2Path; a.Ch3Path};
    end
        
    % load & crop
    LoadCrop (fList, a.cropCoords.*Conditions.ImScale, Conditions); %convert coords from image to nm
    clear cropPos
    
    % ------------------------------------ Density Filtering
    if a.denistyFilter
        display ('Density Filtering Dataset');
        Conditions.distFilt = 200; 
        Conditions.fMode = 1; 
        if isnan (a.densitySD1)
            Conditions.Ch1Filt = [];
        else
            Conditions.Ch1Filt = a.densitySD1;
        end
        if isnan (a.densitySD2)
            Conditions.Ch2Filt = [];
        else
            Conditions.Ch2Filt = a.densitySD2;
        end
        if isnan (a.densitySD3)
            Conditions.Ch3Filt = [];
        else
            Conditions.Ch3Filt = a.densitySD3;
        end
        Conditions.saveTIFF = 0;
        densityFilter ([Conditions.fName ' - cropped.mat'], Conditions);
        Conditions.fName = ['Filtered - ' Conditions.fName];
        
    end
    
    % ------------------------------------ SAA
    if a.SAAcheck
        display (' ');
        
        if a.SAAuserCDC
            Conditions.SD = 0; %0 causes cutoff value to be used
            Conditions.Cutoff = str2double(a.SAAuserCDCvalue);
        else
            Conditions.SD = str2double(a.CDC_SD); %calculate CDC
            Conditions.Cutoff = str2double(a.SAAcutoff);
        end
        Conditions.BinMax = ceil(str2double(a.SAAdist)); 
        Conditions.Iterations = ceil(str2double(a.SAArand));
        Conditions.AnalyzedFraction = str2double(a.SAAaF);
        Conditions.polyRand = uint8 (a.polyRand);
        if Conditions.AnalyzedFraction > 1
            Conditions.AnalyzedFraction = 1;
        elseif Conditions.AnalyzedFraction <= 0
            Conditions.AnalyzedFraction = 0.1;
        end
        Conditions.gFilter = 3;
        Conditions.ContourDensity = 25;

        if a.nCh == 1 && a.SAAcheck
            display ('Insufficient channels for SAA analysis, skipping to next assay');
        elseif a.nCh == 2 && a.SAAcheck
            MIiSRdata.SAAdata = SAA2col([Conditions.fName ' - cropped.mat'], 0, Conditions);
        elseif a.nCh == 3 && a.SAAcheck
            MIiSRdata.SAAdata = SAA3col([Conditions.fName ' - cropped.mat'], 0, Conditions);
        end
    
        pause (0.5); %provide time for graphs to close
        
    end %end SAA
        
    %-------------------------------------- Spatial Stats
    if a.RDFcheck
        display (' ');
        display ('Performing Spatial Statistics');
        Conditions.Bound = str2double(a.spatialDist);
        Conditions.spatialGraph = 1;
        if a.spatialCh1 > a.nCh
            display ('Warning: Channel number selected for the Primary channel is greater than the number of channels in the image.');
            display (['Setting Primary channel to channel ' num2str(a.nCh), '.']);
            a.spatialCh1 = a.nCh;
        end
        if a.spatialCh2 > a.nCh
            display ('Warning: Channel number selected for the Secondary channel is greater than the number of channels in the image.');
            display (['Setting Secondary channel to channel ' num2str(a.nCh), '.']);
            a.spatialCh2 = a.nCh;
        end
        
        %Process
        MIiSRdata.spatialData = spatialStats([Conditions.fName ' - cropped.mat'], a.spatialCh1, a.spatialCh2, 1, Conditions);
        pause (0.5); %provide time for graphs to close
    end %end spatial stats
    
    %-------------------------------------- DBSCAN
    if a.DBSCANcheck 
        display (' ');
        display (['Performing DBSCAN Segmentation of Channel ' num2str(a.DBSCANch) '.']);
        
        if a.DBSCANch > a.nCh
            display (['Selected channel does not exist in image, processing channel ' num2str(a.nCh) ' instead.']);
            a.DBSCANch = a.nCh;
        end
        
        %collect variables
        dK = ceil(str2double(a.DBSCANk));
        if dK < 2
            display('Warning: Minimum cluster size is too small; setting to 2.');
            dK = 2;
        end
        dE = str2double(a.DBSCANe);
        if isnan(dE)
            dE = [];
        end
        
        %run DBSCAN
    try   
        MIiSRdata.DBSCAN = DBSCAN_v([Conditions.fName ' - cropped.mat'], a.DBSCANch, dK, dE, 1, Conditions.fName);
        pause (0.5); %provide time for graphs to close
    catch
        disp('cannot DBSCAN, likely not enough memory!');
        MIiSRdata.DBSCAN = [];
    end
    end %end DBSCAN
    
    %-------------------------------------- OPTICS
    if a.OPTICScheck
        display (' ');
        display (['Performing OPTICS Segmentation of Channel ' num2str(a.OPTICSch) '.']);
        
               
        if a.OPTICSch > a.nCh
            display (['Selected channel does not exist in image, processing channel ' num2str(a.nCh) ' instead.']);
            a.OPTICSch = a.nCh;
        end
        
         %collect variables
        dK = ceil(str2double(a.OPTICSk));
        if dK < 2
            display('Warning: Minimum cluster size is too small; setting to 2.');
            dK = 2;
        end
        
        %run OPTICS
        MIiSRdata.OPTICS = OPTICS([Conditions.fName ' - cropped.mat'], a.OPTICSch, 2, dK);
        
        %run hieartichal segrigation, if selected
        if a.OPTICShier
            display (['   ...Segmenting OPTICS Plot from Channel ' num2str(a.OPTICSch) '.']);
            minSplit = str2double(a.hierRatio);
            
            if minSplit <= 0 || minSplit >= 1.0
                disply ('Warning: Peak Ratio must be between 0 and 1.0. Setting to defulat value (0.75).');
                minSplit = 0.75;
            end
            
            %run heirOPTICS
            MIiSRdata.OPTICSsegment = hierOPTICS (MIiSRdata.OPTICS, dK, minSplit);
        end
        
    end % end OPTICS
       
    
    %-----------------------------SaveData & Cleanup
    save (['MIiSR - ' Conditions.fName '.mat'], 'Conditions', 'MIiSRdata'); 
    
    %prepare for next loop
    clear SAAdata Conditions fList a MIiSRdata dK dE

    %update display and queuelist
    display (' ');
    qList(1) = [];
    set(handles.queueList, 'String', qList);
    pause (0.5); %provide time for GUI to update

    
end %end main processing loop


%-------------------------Final Cleanup
cd (handles.pwd);

%turn off parallel if not started by user
if isempty(b)
    stopPool_Callback(hObject, [], handles);
end

set(handles.startButton, 'Enable', 'inactive');
set(handles.startButton, 'String', 'Done!');
set(handles.startButton, 'ForegroundColor', [0,1,0]);

display('=============================================================');
display('|                                                           |');
display('|           Analysis Complete                               |');
display('|                                                           |');
display('=============================================================');
%beep;
pause (0.5);
close(gcf);

% --- Executes on button press in delROI.
function delROI_Callback(hObject, eventdata, handles)
% hObject    handle to delROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield (handles, 'h')
    delete (handles.h);
    guidata(hObject, handles);
end


% --- Executes on selection change in spatialCh1.
function spatialCh1_Callback(hObject, eventdata, handles)
% hObject    handle to spatialCh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns spatialCh1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spatialCh1


% --- Executes during object creation, after setting all properties.
function spatialCh1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spatialCh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in spatialCh2.
function spatialCh2_Callback(hObject, eventdata, handles)
% hObject    handle to spatialCh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns spatialCh2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spatialCh2


% --- Executes during object creation, after setting all properties.
function spatialCh2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spatialCh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearImage.
function clearImage_Callback(hObject, eventdata, handles)
% hObject    handle to clearImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%reset GUI
set(handles.newROI, 'Enable', 'off');
set(handles.delROI, 'Enable', 'off');
set(handles.clearImage, 'Enable', 'off');
set(handles.Ch1Path, 'String', []);
set(handles.Ch1Name, 'String', 'Name for Channel 1');
set(handles.addCh2, 'Enable', 'off');
set(handles.Ch2Name, 'String', 'Name for Channel 2');
set(handles.Ch2Path, 'String', []);
set(handles.addCh3, 'Enable', 'off');
set(handles.Ch3Path, 'String', []);
set(handles.Ch3Name, 'String', 'Name for Channel 3');
set(handles.addQueue, 'Enable', 'off');

if isfield (handles, 'h');
    delete (handles.h); %clear ROI
end
axes(handles.imDisplay);
imshow([]); %clear image



function hierRatio_Callback(hObject, eventdata, handles)
% hObject    handle to hierRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hierRatio as text
%        str2double(get(hObject,'String')) returns contents of hierRatio as a double


% --- Executes during object creation, after setting all properties.
function hierRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hierRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in denistyFilter.
function denistyFilter_Callback(hObject, eventdata, handles)
% hObject    handle to denistyFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of denistyFilter
tVal = get(handles.denistyFilter, 'Value');
if tVal
    set (handles.densitySD1, 'Enable', 'on');
    set (handles.densitySD2, 'Enable', 'on');
    set (handles.densitySD3, 'Enable', 'on');
    set (handles.densityDist, 'Enable', 'on');
else
    set (handles.densitySD1, 'Enable', 'off');
    set (handles.densitySD2, 'Enable', 'off');
    set (handles.densitySD3, 'Enable', 'off');
    set (handles.densityDist, 'Enable', 'off');
end


function densitySD2_Callback(hObject, eventdata, handles)
% hObject    handle to densitySD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of densitySD2 as text
%        str2double(get(hObject,'String')) returns contents of densitySD2 as a double


% --- Executes during object creation, after setting all properties.
function densitySD2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to densitySD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function densityDist_Callback(hObject, eventdata, handles)
% hObject    handle to densityDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of densityDist as text
%        str2double(get(hObject,'String')) returns contents of densityDist as a double


% --- Executes during object creation, after setting all properties.
function densityDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to densityDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function densitySD1_Callback(hObject, eventdata, handles)
% hObject    handle to densitySD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of densitySD1 as text
%        str2double(get(hObject,'String')) returns contents of densitySD1 as a double


% --- Executes during object creation, after setting all properties.
function densitySD1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to densitySD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function densitySD3_Callback(hObject, eventdata, handles)
% hObject    handle to densitySD3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of densitySD3 as text
%        str2double(get(hObject,'String')) returns contents of densitySD3 as a double


% --- Executes during object creation, after setting all properties.
function densitySD3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to densitySD3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SAAaF_Callback(hObject, eventdata, handles)
% hObject    handle to SAAaF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SAAaF as text
%        str2double(get(hObject,'String')) returns contents of SAAaF as a double


% --- Executes during object creation, after setting all properties.
function SAAaF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SAAaF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in polyRand.
function polyRand_Callback(hObject, eventdata, handles)
% hObject    handle to polyRand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of polyRand
