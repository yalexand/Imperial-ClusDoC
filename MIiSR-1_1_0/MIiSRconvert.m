function varargout = MIiSRconvert(varargin)
% MIISRCONVERT MATLAB code for MIiSRconvert.fig
%      MIISRCONVERT, by itself, creates a new MIISRCONVERT or raises the existing
%      singleton*.
%
%      H = MIISRCONVERT returns the handle to a new MIISRCONVERT or the handle to
%      the existing singleton*.
%
%      MIISRCONVERT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MIISRCONVERT.M with the given input arguments.
%
%      MIISRCONVERT('Property','Value',...) creates a new MIISRCONVERT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MIiSRconvert_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MIiSRconvert_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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
% Edit the above text to modify the response to help MIiSRconvert

% Last Modified by GUIDE v2.5 09-Jul-2015 10:34:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MIiSRconvert_OpeningFcn, ...
                   'gui_OutputFcn',  @MIiSRconvert_OutputFcn, ...
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


% --- Executes just before MIiSRconvert is made visible.
function MIiSRconvert_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MIiSRconvert (see VARARGIN)

% Choose default command line output for MIiSRconvert
handles.output = hObject;
handles.fileList = []; %structured: path/file, file type, photon filter, image scale, intensity scale

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MIiSRconvert wait for user response (see UIRESUME)
% uiwait(handles.MIiSRconv);

% --- Outputs from this function are returned to the command line.
function varargout = MIiSRconvert_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function binSize_Callback(hObject, eventdata, handles)
% hObject    handle to binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binSize as text
%        str2double(get(hObject,'String')) returns contents of binSize as a double

binSize = ceil(str2num(get(handles.binSize, 'String')));
if binSize < 1
    binSize = 1;
end
binSize = num2str(binSize);
set (handles.binSize, 'String', num2str(binSize));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function binSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iFilterPhoton_Callback(hObject, eventdata, handles)
% hObject    handle to iFilterPhoton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iFilterPhoton as text
%        str2double(get(hObject,'String')) returns contents of iFilterPhoton as a double
iFilter = str2num(get(handles.iFilterPhoton, 'String'));
if iFilter < 0;
    iFilter = 0;
end

handles.iFilterPhoton = ceil(iFilter); %ensure iFilter is an intiger
set (handles.iFilterPrecision, 'String', num2str(200/sqrt(handles.iFilterPhoton), '%f'));

% --- Executes during object creation, after setting all properties.
function iFilterPhoton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iFilterPhoton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iFilterPrecision_Callback(hObject, eventdata, handles)
% hObject    handle to iFilterPrecision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iFilterPrecision as text
%        str2double(get(hObject,'String')) returns contents of iFilterPrecision as a double
iFilter = str2num(get(handles.iFilterPrecision, 'String'));
if iFilter < 0;
    iFilter = 0;
end
handles.iFilterPrecision = iFilter;
set (handles.iFilterPhoton, 'String', num2str((200/iFilter)^2));


% --- Executes during object creation, after setting all properties.
function iFilterPrecision_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iFilterPrecision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Scale_Callback(hObject, eventdata, handles)
% hObject    handle to Scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Scale as text
%        str2double(get(hObject,'String')) returns contents of Scale as a double


% --- Executes during object creation, after setting all properties.
function Scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function iScale_Callback(hObject, eventdata, handles)
% hObject    handle to iScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iScale as text
%        str2double(get(hObject,'String')) returns contents of iScale as a double


% --- Executes during object creation, after setting all properties.
function iScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in fileType.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to fileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fileType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fileType


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in fileType.
function fileType_Callback(hObject, eventdata, handles)
% hObject    handle to fileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fileType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fileType


% --- Executes during object creation, after setting all properties.
function fileType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addFile.
function addFile_Callback(hObject, eventdata, handles)
% hObject    handle to addFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.StartRun, 'Enable', 'off');
set(handles.addSubFolder, 'Enable', 'off');
set(handles.addFile, 'Enable', 'off');
set(handles.addFolder, 'Enable', 'off');

fType = get(handles.fileType, 'Value'); %get selected file type
if fType == 1
    fExt = '.xls';
elseif fType == 2
    fExt = '.ascii';
else
    fExt = '.txt';
end

fGet = ['*' fExt];

[tmpFile, tmpPath] = uigetfile(fGet, ['Select A Single ' fExt 'File']); % get file
qList = get(handles.queueList, 'String'); %get current display queue
fileList = handles.fileList; %get current processing queue

if length(tmpFile)>1
    tmpFile = [tmpPath, tmpFile];
    tmpFile = cellstr(tmpFile);
    qList = [qList; tmpFile];
    set(handles.queueList, 'String', qList);
    tFilter = str2double(get(handles.iFilterPhoton, 'String'));
    tScale = str2double(get(handles.Scale, 'String'));
    tiScale = str2double(get(handles.iScale, 'String'));
    
    tFiles = [tmpFile, num2cell(fType), num2cell(tFilter), num2cell(tScale), num2cell(tiScale)];
    handles.fileList = [fileList; tFiles];
    guidata(hObject, handles);
else
    display ('No files found in the selected directory');
end

set(handles.StartRun, 'Enable', 'on');
set(handles.addSubFolder, 'Enable', 'on');
set(handles.addFile, 'Enable', 'on');
set(handles.addFolder, 'Enable', 'on');

% --- Executes on button press in addFolder.
function addFolder_Callback(hObject, eventdata, handles)
% hObject    handle to addFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.StartRun, 'Enable', 'off');
set(handles.addSubFolder, 'Enable', 'off');
set(handles.addFile, 'Enable', 'off');
set(handles.addFolder, 'Enable', 'off');

fType = get(handles.fileType, 'Value'); %get selected file type
if fType == 1
    fExt = '.xls';
elseif fType == 2
    fExt = '.ascii';
else
    fExt = '.txt';
end


tmpDir = uigetdir('', ['Select All ' fExt 'Files in a Single Directory']); % get directory
qList = get(handles.queueList, 'String'); %get current display queue
fileList = handles.fileList; %get current processing queue

if length(tmpDir) > 1
    tmpDir = cellstr(tmpDir);
	tFiles = getFiles(char(tmpDir), fExt);
    qList = [qList; tFiles];
    set(handles.queueList, 'String', qList);
    tFilter = str2double(get(handles.iFilterPhoton, 'String'));
    tScale = str2double(get(handles.Scale, 'String'));
    tiScale = str2double(get(handles.iScale, 'String'));
    
    for ii=1:length(tFiles)
        tFiles(ii,2) = num2cell(fType);
        tFiles(ii,3) = num2cell(tFilter);
        tFiles(ii,4) = num2cell(tScale);
        tFiles(ii,5) = num2cell(tiScale);
    end
    handles.fileList = [fileList; tFiles];
    guidata(hObject, handles);
else
    display ('No files found in the selected directory');
end

set(handles.StartRun, 'Enable', 'on');
set(handles.addSubFolder, 'Enable', 'on');
set(handles.addFile, 'Enable', 'on');
set(handles.addFolder, 'Enable', 'on');

% --- Executes on button press in addSubFolder.
function addSubFolder_Callback(hObject, eventdata, handles)
% hObject    handle to addSubFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.StartRun, 'Enable', 'off');
set(handles.addSubFolder, 'Enable', 'off');
set(handles.addFile, 'Enable', 'off');
set(handles.addFolder, 'Enable', 'off');

fType = get(handles.fileType, 'Value'); %get selected file type
if fType == 1
    fExt = '.xls';
elseif fType == 2
    fExt = '.ascii';
else
    fExt = '.txt';
end

tmpDir = uigetdir('', ['Select All ' fExt 'Files in All Subfolders']); % get directories
qList = get(handles.queueList, 'String'); %get current display queue
fileList = handles.fileList; %get current processing queue
tFilter = str2double(get(handles.iFilterPhoton, 'String'));
tScale = str2double(get(handles.Scale, 'String'));
tiScale = str2double(get(handles.iScale, 'String'));

tList = [];
if length(tmpDir) > 1
    folderList = foldersRecurs(tmpDir, fExt);
    for ii=1:size(folderList,1)
        tFiles = getFiles(char(folderList(ii,:)), fExt);
        tList = [tList;tFiles];
    end
    
    qList = [qList; tList];
    set(handles.queueList, 'String', qList);
    
    for ii=1:length(tList)
        tList(ii,2) = num2cell(fType);
        tList(ii,3) = num2cell(tFilter);
        tList(ii,4) = num2cell(tScale);
        tList(ii,5) = num2cell(tiScale);
    end
    handles.fileList = [fileList; tList];
    guidata(hObject, handles);
else
    display ('No files found in the selected directory');
end

set(handles.StartRun, 'Enable', 'on');
set(handles.addSubFolder, 'Enable', 'on');
set(handles.addFile, 'Enable', 'on');
set(handles.addFolder, 'Enable', 'on');


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
HelpFile = which('MIiSRconvert.m');
s=strfind(HelpFile,'MIiSRconvert.m');
HelpFile(s:end)=[];
HelpFile = [HelpFile 'help' filesep 'MIiSRconvertHelp.html'];
web (HelpFile);


% --- Executes on button press in delEntry.
function delEntry_Callback(hObject, eventdata, handles)
% hObject    handle to delEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selLine = get(handles.queueList, 'Value');
qList = get(handles.queueList, 'String'); %get current display queue
fileList = handles.fileList; %get current processing queue
qList(selLine) = [];
fileList(selLine,:) = [];
set(handles.queueList, 'String', qList);
handles.fileList = fileList;
guidata(hObject, handles);


function queueList_CreateFcn(hObject, eventdata, handles)

% --- Executes during object deletion, before destroying properties.
function queueList_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to queueList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in StartRun.
function StartRun_Callback(hObject, eventdata, handles)
% hObject    handle to StartRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set GUI to inactive
guidata(hObject, handles);

set(handles.StartRun, 'Enable', 'off');
set(handles.addSubFolder, 'Enable', 'off');
set(handles.addFile, 'Enable', 'off');
set(handles.addFolder, 'Enable', 'off');
set(handles.binSize, 'Enable', 'off');
set(handles.iFilterPrecision, 'Enable', 'off');
set(handles.iFilterPhoton, 'Enable', 'off');
set(handles.Scale, 'Enable', 'off');
set(handles.helpButton, 'Enable', 'off');
set(handles.queueList, 'Enable', 'off');
set(handles.delEntry, 'Enable', 'off');
set(handles.fileType, 'Enable', 'off');
set(handles.iScale, 'Enable', 'off');

%collect variables
pause (0.5);
binSize = str2double(get(handles.binSize,'String'));
fList = handles.fileList;

cDir = pwd; %store current directory; return at end

%Process each file
for ii=1:size(fList,1)
    %mark file currently being processed
    qList = get(handles.queueList, 'String'); %get current display queue
    tName = char(qList(1));
    tName = cellstr(['Processing: ' tName]);
    qList(1) = tName;
    set(handles.queueList, 'String', qList);
    
    %find & move to directory
    aDir = char(fList(ii,1));
    aFile = aDir;
    s = strfind(aDir, filesep);
    s = s(end);
    aDir(s:end) = [];
    aFile(1:s) = [];
    cd (aDir);
    pause(1.0);
    
    %run conversion
    Conditions.dataType = cell2mat(fList(ii,2));
    Conditions.iFilter = cell2mat(fList(ii,3));
    Conditions.cropRange = [];
    Conditions.Scale = cell2mat(fList(ii,4));
    Conditions.iScale = cell2mat(fList(ii,5));
    fileConv (aFile, binSize, Conditions);
    
    qList(1) = [];
    set(handles.queueList, 'String', qList);
end

cd (cDir);

set(handles.StartRun, 'Enable', 'on');
set(handles.addSubFolder, 'Enable', 'on');
set(handles.addFile, 'Enable', 'on');
set(handles.addFolder, 'Enable', 'on');
set(handles.binSize, 'Enable', 'on');
set(handles.iFilterPrecision, 'Enable', 'on');
set(handles.iFilterPhoton, 'Enable', 'on');
set(handles.Scale, 'Enable', 'on');
set(handles.helpButton, 'Enable', 'on');
set(handles.queueList, 'Enable', 'on');
set(handles.delEntry, 'Enable', 'on');
set(handles.fileType, 'Enable', 'on');
set(handles.iScale, 'Enable', 'on');

display ('Conversion Complete');
beep;
close(gcf);

% --- Executes on selection change in queueList.
function queueList_Callback(hObject, eventdata, handles)
% hObject    handle to queueList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns queueList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from queueList

%% ===== Additional Functions ===
function files = getFiles(Folder, Extention)
% returns a list of the files and the folders in the location specified by the
% input argument 'Folder'
a = strfind(Folder, filesep);

files=[];
f=dir(Folder);

% set counters
x=1;

while x <= length(f)
    if f(x).isdir
        f(x)=[];
    elseif strcmp(f(x).name,'.') || strcmp(f(x).name,'..')
        f(x)=[];
    else        
        files{x}=f(x).name;
        x=x+1;
    end
end

%remove mac's idiotic additional files
if ~isempty(files)
    fNum = 0;
    for i=1:length(files)
        tName = char(files(i));
        if tName(1) ~= '.'
            fNum = fNum + 1;
            tList(fNum) = files(i);
        end
    end
    if fNum == 0
      files = [];
    else
      files = tList;
    end
end


%reduces list to files present in folder with desired extention
if ~isempty(Extention) && ~isempty(files)
    file_count = size(files,2);
    nfiles = 0;
    for k=1:file_count
        %generate list of files to process
        CurFile = char(files(k));
        s=strfind(CurFile,Extention);
        if ~(isempty(s))
            nfiles = nfiles+1;
            fileList(nfiles,:) = cellstr([Folder, filesep, CurFile]);
        end %endif
    end %end for (k)
    if nfiles == 0
        fileList = char([]);
    end
    clear files;
    files(:,:) = fileList(:,:);
end


function folderStruct = foldersRecurs (homeFolder, fileType)
%
% Function which returns a list of all sub-folders within a home folder
% containing files of the expected file type
%
contents = dir(homeFolder);
directories = find([contents.isdir]);
nFolders = 0;
tmpStruct = num2cell([]);

% For loop will be skipped when directory contains no sub-directories
for i_dir = directories
    sub_directory  = contents(i_dir).name;
    full_directory = fullfile(homeFolder, sub_directory);

    % ignore '.' and '..'
    if (strcmp(sub_directory, '.') || strcmp(sub_directory, '..'))
        continue;
    end
    
    nFolders = nFolders + 1;
    tmpStruct(nFolders ,1) = cellstr(full_directory);
    
    % Recurse down
    tmpStruct2 = detectFoldersRecursive (full_directory, fileType);
    tLen = length(tmpStruct2);
    if (tLen > 0)
       nFolders = nFolders+1;
       tmpStruct(nFolders:nFolders+tLen-1, 1) = tmpStruct2;
       nFolders = nFolders + tLen - 1;
    end
    clear tmpStruct2 tLen;
end

% remove folders lacking filetype
fLen = length(tmpStruct);
fPos = 1;
fStruct = num2cell([]);
for i=1:fLen
    foldersX = getFiles(char(tmpStruct(i)), fileType);
    if ~isempty(foldersX)
        fStruct(fPos,1) = tmpStruct(i);
        fPos = fPos + 1;
    end
end
folderStruct = fStruct;
