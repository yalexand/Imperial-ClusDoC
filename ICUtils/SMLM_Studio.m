function varargout = SMLM_Studio(varargin)
% SMLM_STUDIO MATLAB code for SMLM_Studio.fig
%      SMLM_STUDIO, by itself, creates a new SMLM_STUDIO or raises the existing
%      singleton*.
%
%      H = SMLM_STUDIO returns the handle to a new SMLM_STUDIO or the handle to
%      the existing singleton*.
%
%      SMLM_STUDIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMLM_STUDIO.M with the given input arguments.
%
%      SMLM_STUDIO('Property','Value',...) creates a new SMLM_STUDIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SMLM_Studio_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SMLM_Studio_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SMLM_Studio

% Last Modified by GUIDE v2.5 20-Jan-2022 18:05:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SMLM_Studio_OpeningFcn, ...
                   'gui_OutputFcn',  @SMLM_Studio_OutputFcn, ...
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


% --- Executes just before SMLM_Studio is made visible.
function SMLM_Studio_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SMLM_Studio (see VARARGIN)

% Choose default command line output for SMLM_Studio
handles.output = hObject;

handles.figureName = get(handles.figure1,'Name');

handles.data = [];
handles.filenames = [];

handles.param_names = {'cl.DoC','cl.Area','cl.Circularity','cl.Nloc','cl.Density','roi.Ripley','roi.xRDF','roi.xRipley','roi.Ripley(MAX)','roi.xRDF(MAX)','roi.xRipley(MAX)','roi.SAA','roi.SAA(ratio)' 'roi.ClusterDensity'};

set(handles.axes1,'XTick',[]);
set(handles.axes1,'YTick',[]);
set(handles.corr_plot,'XTick',[]);
set(handles.corr_plot,'YTick',[]);
set(handles.corr_plot2,'XTick',[]);
set(handles.corr_plot2,'YTick',[]);

handles.statistics_names = {'p-value: KS','p-value: t-test','p-value: Wilcoxon','Cohen"s d','|median diff|'};

set(handles.bar_log,'Enable','off');
set(handles.bar2_log,'Enable','off');

condition = {...
    'C1', ...
    'C2', ...
    'C3', ...    
    'C4', ...
    'C5'};
assert(length(unique(condition))==length(condition));
%
n_conditions = length(condition);
n_plates = 2;
n_wells = 7;
n_fovs = 3;
n_channels = 2;
n_max_objects = 6;          

% red chooser
   data = [num2cell(true(n_conditions,1)) condition'];
   set(handles.Condition_table, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
         
   data = [num2cell(true(n_plates,1)) num2cell((1:n_plates))'];
   set(handles.Plate_table, 'Data', data, ... 
         'ColumnEditable',[true,false]);  

   data = [num2cell(true(n_wells,1)) num2cell((1:n_wells))'];
   set(handles.Well_table, 'Data', data, ... 
         'ColumnEditable',[true,false]);  

   data = [num2cell(true(n_fovs,1)) num2cell((1:n_fovs))'];
   set(handles.FOV_table, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
     
   data = [num2cell(true(n_channels,1)) num2cell((1:n_channels))'];
   set(handles.Channel_table, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
     
% blue chooser
   data = [num2cell(true(n_conditions,1)) condition'];
   set(handles.Condition_table2, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
     
   data = [num2cell(true(n_plates,1)) num2cell((1:n_plates))'];
   set(handles.Plate_table2, 'Data', data, ... 
         'ColumnEditable',[true,false]);  

   data = [num2cell(true(n_wells,1)) num2cell((1:n_wells))'];
   set(handles.Well_table2, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
   data = [num2cell(true(n_fovs,1)) num2cell((1:n_fovs))'];
   set(handles.FOV_table2, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
     
   data = [num2cell(true(n_channels,1)) num2cell((1:n_channels))'];
   set(handles.Channel_table2, 'Data', data, ... 
         'ColumnEditable',[true,false]);  

% 2d histo chooser
   data = [num2cell(true(n_conditions,1)) condition'];
   set(handles.Condition_table3, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
     
   data = [num2cell(true(n_plates,1)) num2cell((1:n_plates))'];
   set(handles.Plate_table3, 'Data', data, ... 
         'ColumnEditable',[true,false]);  

   data = [num2cell(true(n_wells,1)) num2cell((1:n_wells))'];
   set(handles.Well_table3, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
   data = [num2cell(true(n_fovs,1)) num2cell((1:n_fovs))'];
   set(handles.FOV_table3, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
     
   data = [num2cell(true(n_channels,1)) num2cell((1:n_channels))'];
   set(handles.Channel_table3, 'Data', data, ... 
         'ColumnEditable',[true,false]);  

    set(handles.Condition_table,'CellEditCallback',@table_check_callback);
        set(handles.Condition_table2,'CellEditCallback',@table_check_callback);
            set(handles.Condition_table3,'CellEditCallback',@table_check_callback);
     
    set(handles.Plate_table,'CellEditCallback',@table_check_callback);
        set(handles.Plate_table2,'CellEditCallback',@table_check_callback);
            set(handles.Plate_table3,'CellEditCallback',@table_check_callback);

    set(handles.Well_table,'CellEditCallback',@table_check_callback);
        set(handles.Well_table2,'CellEditCallback',@table_check_callback);
            set(handles.Well_table3,'CellEditCallback',@table_check_callback);

    set(handles.FOV_table,'CellEditCallback',@table_check_callback);
        set(handles.FOV_table2,'CellEditCallback',@table_check_callback);
            set(handles.FOV_table3,'CellEditCallback',@table_check_callback);

    set(handles.Channel_table,'CellEditCallback',@table_check_callback);
        set(handles.Channel_table2,'CellEditCallback',@table_check_callback);
            set(handles.Channel_table3,'CellEditCallback',@table_check_callback);

    set(handles.Object_table,'CellEditCallback',@table_check_callback);
        set(handles.Object_table2,'CellEditCallback',@table_check_callback);
            set(handles.Object_table3,'CellEditCallback',@table_check_callback);
           
%%%%%%%
% 2d histo chooser
   data = [num2cell(true(n_conditions,1)) condition'];
   set(handles.Condition_table4, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
     
   data = [num2cell(true(n_plates,1)) num2cell((1:n_plates))'];
   set(handles.Plate_table4, 'Data', data, ... 
         'ColumnEditable',[true,false]);  

   data = [num2cell(true(n_wells,1)) num2cell((1:n_wells))'];
   set(handles.Well_table4, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
   data = [num2cell(true(n_fovs,1)) num2cell((1:n_fovs))'];
   set(handles.FOV_table4, 'Data', data, ... 
         'ColumnEditable',[true,false]);  
     
   data = [num2cell(true(n_channels,1)) num2cell((1:n_channels))'];
   set(handles.Channel_table4, 'Data', data, ... 
         'ColumnEditable',[true,false]);  

            set(handles.Condition_table4,'CellEditCallback',@table_check_callback);     
            set(handles.Plate_table4,'CellEditCallback',@table_check_callback);
            set(handles.Well_table4,'CellEditCallback',@table_check_callback);
            set(handles.FOV_table4,'CellEditCallback',@table_check_callback);
            set(handles.Channel_table4,'CellEditCallback',@table_check_callback);
            set(handles.Object_table4,'CellEditCallback',@table_check_callback);
%%%%%%%
            
            
handles.data = cell( ...
            n_plates, ...
            n_conditions, ...
            n_wells, ...
            n_fovs, ...
            n_channels, ...
            n_max_objects);  
        
handles.mask = []; % no data 
handles.tot_data = [];
%
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SMLM_Studio wait for user response (see UIRESUME)
%uiwait(handles.figure1);

function table_check_callback(hObject,callbackdata)

     handles = guidata(hObject);

     if isempty(handles.tot_data)
         disp('not doing any work, data were not loaded');
         return;
     end
    
     Tag = get(hObject,'Tag');
     
     postfix = {'_table','_table2','_table3','_table4'};
     
     if contains(Tag,postfix{2})
         k = 2;
     elseif contains(Tag,postfix{3})
         k = 3;
     elseif contains(Tag,postfix{4})
         k = 4;         
     else
         k = 1;
     end
        P = eval(['handles.Plate' postfix{k}]);        
        C = eval(['handles.Condition' postfix{k}]);
        W = eval(['handles.Well' postfix{k}]);
        F = eval(['handles.FOV' postfix{k}]);  
        O = eval(['handles.Object' postfix{k}]);          
        c = eval(['handles.Channel' postfix{k}]);
        %
        d = get(P,'Data');
        iP = d(:,1);
            d = get(C,'Data');        
            iC = d(:,1);
                d = get(W,'Data');
                iW = d(:,1);
                    d = get(F,'Data');
                    iF = d(:,1);
                        d = get(O,'Data');
                        iO = d(:,1);                    
                            d = get(c,'Data');
                            ic = d(:,1);
                            
     %
    cur_P = handles.Plate(cell2mat(iP));
    cur_C = handles.Condition(cell2mat(iC));    
    % !!! cur_W = handles.Well(cell2mat(t{3}));
        uW = (handles.unique_wells)';
        cur_W = uW(cell2mat(iW));
    %
    cur_F = handles.FOVs(cell2mat(iF));
    cur_O = handles.unique_objects(cell2mat(iO));
    cur_c = handles.Channel(cell2mat(ic));
                                 
    N = numel(handles.tot_data);
    
    for m=1:N
        [P,C,W,F,O,~,c] = get_attributes(handles.tot_data{m}.token);
        %
        handles.mask(m,k) = 1;               
        %
        if 0==sum(ismember(cur_P,{P}))
            handles.mask(m,k) = 0;
            continue;
        end
        if 0==sum(ismember(cur_C,{C}))
            handles.mask(m,k) = 0;
            continue;
        end
        if 0==sum(ismember(cur_W,{W}))
            handles.mask(m,k) = 0;
            continue;
        end        
        if 0==sum(ismember(cur_c,{c}))
            handles.mask(m,k) = 0;
            continue;
        end  
        %
        fov_token = [':P:' P ':C:' C  ':W:' W  ':F:' F];
        if 0==sum(ismember(cur_F,fov_token))
            handles.mask(m,k) = 0;
            continue;
        end  
        obj_token = [fov_token ':O:' O];
        if 0==sum(ismember(cur_O,obj_token))
            handles.mask(m,k) = 0;
            continue;
        end                          
    end
    sum(handles.mask,1)
    %
    guidata(hObject,handles);
    if ismember(k,[1 2])
        show_plot(hObject,handles);
    elseif 3==k
        show_2d_histogram(handles);
    else
        show_2d_histogram2(handles);        
    end
          
% --- Outputs from this function are returned to the command line.
function varargout = SMLM_Studio_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in bar_histo.
function bar_histo_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    logflag = 'on';
else
    logflag = 'off';
end
    set(handles.bar_log,'Enable',logflag);
show_2d_histogram(handles);


% --- Executes on button press in bar_log.
function bar_log_Callback(hObject, eventdata, handles)
show_2d_histogram(handles);

function minmaxlimits = define_minmax(handles)
minmaxlimits = zeros(numel(handles.param_names),2);


function Q_X_Callback(hObject, eventdata, handles)
show_2d_histogram(handles); 

% --- Executes during object creation, after setting all properties.
function Q_X_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Q_Y_Callback(hObject, eventdata, handles)
show_2d_histogram(handles); 

% --- Executes during object creation, after setting all properties.
function Q_Y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Q1_Callback(hObject, eventdata, handles)
show_plot(hObject,handles); 

% --- Executes during object creation, after setting all properties.
function Q1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % % --------------------------------------------------------------------
function load_data_Callback(hObject, eventdata, handles)
%
[filename,pathname] = uigetfile('*.mat','Select data file',pwd);                
if isempty(filename), return, end       
if isnumeric(filename) && 0==filename, return, end

handles.fullfilename = [pathname filesep filename];
load(handles.fullfilename);

well_to_condition_map = SMLMdata.well_to_condition_map;
well_to_condition_index_map = SMLMdata.well_to_condition_index_map;

           [n_plates, ...
            n_conditions, ...
            n_wells, ...
            n_fovs, ...
            n_channels, ...
            n_max_objects] = size(SMLMdata.data);
tic

%cnt_objects = 0;

handles.tot_data = cell(1,1000000); % million ROIs

cnt_rois = 0;

for plate = 1:n_plates
    for w=1:n_wells        
        well_is_empty = true;
        try
            well_is_empty = isempty(well_to_condition_map(num2str(w)));
        catch
            well_is_empty = isempty(well_to_condition_map(SMLMdata.Well{w}));
        end
        %
        fov = w; % this is lkely wrong, but held for first data reduction cases
        %
        for k = 1:n_max_objects
            for channel = 1:n_channels
                %
                try
                    condition_index = well_to_condition_index_map(num2str(w));
                catch
                    condition_index = well_to_condition_index_map(SMLMdata.Well{w});
                end
                %
                cur_data = SMLMdata.data{plate,condition_index,w,fov,channel,k};
                if isempty(cur_data), continue, end
                nrois = numel(cur_data);
                    for r=1:nrois
                        ass = cur_data{r};
                        cnt_rois = cnt_rois+1;
                        handles.tot_data{cnt_rois} = ass;
                    end
            end                                
        end
    end
end
%
handles.tot_data = handles.tot_data(1,1:cnt_rois); % less than million..

% check if there are centroid coordinates

there_are_cluster_centroid_coordinates = true;
        ass = handles.tot_data{1};
        if isfield(ass,'DBSCAN_clusters') && ~isempty(ass.DBSCAN_clusters{1})
                if ~isfield(ass.DBSCAN_clusters{1},'Xc')
                    there_are_cluster_centroid_coordinates = false;                     
                end
        end
%        
Lmax = 300;
if there_are_cluster_centroid_coordinates
    handles.tot_data = get_cluster_density_per_ROI(handles,Lmax);
    guidata(hObject,handles);
end

toc/60

handles.Plate = SMLMdata.Plate;
handles.data = SMLMdata.data;
handles.Well = SMLMdata.Well;
handles.Condition = SMLMdata.Condition;
handles.FOVs = SMLMdata.FOVs;
handles.Object = SMLMdata.Object;
handles.Channel = SMLMdata.Channel;
handles.well_to_condition_map = SMLMdata.well_to_condition_map;
handles.well_to_condition_index_map = SMLMdata.well_to_condition_index_map;
handles.Ripley_distance = SMLMdata.Ripley_distance;

data = [num2cell(true(n_plates,1)) num2cell((1:n_plates))'];
%data(:,2) = (handles.Plate)'; %? error
set(handles.Plate_table, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Plate_table2, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Plate_table3, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Plate_table4, 'Data', data, 'ColumnEditable',[true,false]);
%
data = [num2cell(true(n_conditions,1)) num2cell((1:n_conditions))'];
data(:,2) = (handles.Condition)';
set(handles.Condition_table, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Condition_table2, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Condition_table3, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Condition_table4, 'Data', data, 'ColumnEditable',[true,false]);

data = [num2cell(true(n_channels,1)) num2cell((1:n_channels))'];
data(:,2) = (handles.Channel)';
set(handles.Channel_table, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Channel_table2, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Channel_table3, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Channel_table4, 'Data', data, 'ColumnEditable',[true,false]);

% less trivial
%FOVs
data = [num2cell(true(numel(handles.FOVs),1)) num2cell((1:numel(handles.FOVs)))'];
data(:,2) = handles.FOVs;
set(handles.FOV_table, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.FOV_table2, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.FOV_table3, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.FOV_table4, 'Data', data, 'ColumnEditable',[true,false]);

%Wells
unique_wells = cell(0);
for k = 1:n_wells
    if ~isempty(well_to_condition_map(handles.Well{k}))
        unique_wells = [unique_wells; handles.Well{k}];
    end
end
n_unique_wells = numel(unique_wells);
data = [num2cell(true(n_unique_wells,1)) num2cell((1:n_unique_wells))'];
data(:,2) = unique_wells;
set(handles.Well_table, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Well_table2, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Well_table3, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Well_table4, 'Data', data, 'ColumnEditable',[true,false]);
handles.unique_wells = unique_wells;

%Objects
these_ones = [];
for k=1:numel(handles.tot_data)
    [P,C,W,F,O,~,~] = get_attributes(handles.tot_data{k}.token);
    these_ones = [these_ones; {[':P:' P ':C:' C ':W:' W ':F:' F ':O:' O]}];
end
handles.unique_objects = unique(these_ones);
n_unique_objects = numel(handles.unique_objects);
data = [num2cell(true(n_unique_objects,1)) num2cell((1:n_unique_objects))'];
data(:,2) = handles.unique_objects;
set(handles.Object_table, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Object_table2, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Object_table3, 'Data', data, 'ColumnEditable',[true,false]);
set(handles.Object_table4, 'Data', data, 'ColumnEditable',[true,false]);
%Objects

handles.mask = ones(numel(handles.tot_data),4); % mask for 4 choosers

handles.param_names = cell(0);

there_are_clusters = false;

cluster_params = {'Area','Nb','Circularity','Elongation','MeanDoC'};
for k=1:numel(handles.tot_data)
        ass = handles.tot_data{k};
        if isfield(ass,'DBSCAN_clusters') && ~isempty(ass.DBSCAN_clusters{1})
            for m=1:numel(cluster_params)
                param_name = cluster_params{m};
                if isfield(ass.DBSCAN_clusters{1},param_name)
                    handles.param_names = [handles.param_names; ['cl.' param_name]];
                end
            end
            handles.param_names = [handles.param_names; 'cl.Density'; 'cl.NormDensity']; % to be calculated on the fly
            there_are_clusters = true;
            break;
        end
end

there_are_Ripley_data = false;
for k=1:numel(handles.tot_data)
    ass = handles.tot_data{k};
    if isfield(ass,'RipleyK_val')
        handles.param_names = [handles.param_names; 'roi.Ripley'; 'roi.Ripley(MAX)'];
        there_are_Ripley_data = true;
        break;
    end
end

if there_are_cluster_centroid_coordinates
    handles.param_names  = [handles.param_names; 'roi.ClusterDensity'];
    guidata(hObject,handles);
end

set(handles.Q1,'String',handles.param_names);
set(handles.Q_X,'String',handles.param_names);
set(handles.Q_Y,'String',handles.param_names);
set(handles.Q2_X,'String',handles.param_names);
set(handles.Q2_Y,'String',handles.param_names);
%
% needs normalizing clusters' locaization density 
% by FOV's average localizations density.
handles.FOV_ind = containers.Map(handles.FOVs,1:length(handles.FOVs));
handles.Chan_ind = containers.Map(handles.Channel,1:length(handles.Channel));
%
if there_are_clusters
    N_acc = zeros(length(handles.FOVs),2); % 2 channels
    A_acc = zeros(length(handles.FOVs),2);
        for k=1:numel(handles.tot_data)
                ass = handles.tot_data{k};
                if isfield(ass,'DBSCAN_clusters') && ~isempty(ass.DBSCAN_clusters)
                    [P,C,W,F,~,~,c] = get_attributes(ass.token);
                    index = handles.FOV_ind([':P:' P ':C:' C ':W:' W ':F:' F]);
                    cind = handles.Chan_ind(c);
                    for cl = 1:numel(ass.DBSCAN_clusters)
                        N_acc(index,cind) = N_acc(index,cind) + ass.DBSCAN_clusters{cl}.Nb;
                        A_acc(index,cind) = A_acc(index,cind) + ass.DBSCAN_clusters{cl}.Area;
                    end
                end
        end
handles.FOVs_density = N_acc./A_acc;
end

if there_are_Ripley_data
    step_nm = handles.Ripley_distance(2)-handles.Ripley_distance(1);
    handles.RipleyK_data_MAX = nan(1,numel(handles.tot_data));
            for k=1:numel(handles.tot_data)
                Y = handles.tot_data{k}.RipleyK_val;
                if ~isempty(Y)
                    maxind = find(Y==max(Y));
                    maxind = maxind(1); % :)                    
                    handles.RipleyK_data_MAX(k) = (maxind-1)*step_nm;
                end
            end
end

% set up minmax
handles.minmaxlimits = zeros(numel(handles.param_names),2); 
for k=1:numel(handles.param_names)        
    tic
    set(handles.Q1,'Value',k);
    [s, Q, ~] = select_data(handles,'Q1',1);    
    Q    
    try
    handles.minmaxlimits(k,:) = [min(s) max(s)];
    catch
        disp('error!')
        disp(handles.param_names{k});
    end
    toc
    switch handles.param_names{k}
        case 'cl.Area'
        case 'cl.Nb'
        case 'cl.Circularity'
            handles.minmaxlimits(k,:) = [0 1];            
        case 'cl.Elongation'
            handles.minmaxlimits(k,:) = [0 1];            
        case 'cl.MeanDoC'
            handles.minmaxlimits(k,:) = [-1 1];            
        case 'cl.Density'
        case 'cl.NormDensity'
        case 'roi.Ripley'
        case 'roi.Ripley(MAX)'
        case 'roi.ClusterDensity'
    end
end

%
% Update handles structure
guidata(hObject, handles);

set(handles.figure1, 'Name', [handles.figureName ' : ' handles.fullfilename]);

show_plot(hObject,handles);
show_2d_histogram(handles);
show_2d_histogram2(handles);

%--------------------------------------------------------------------
function [P,C,W,F,O,R,c] = get_attributes(token)
P = extractBetween(token,':P:',':C:');
C = extractBetween(token,':C:',':W:');
W = extractBetween(token,':W:',':F:');
F = extractBetween(token,':F:',':O:');
O = extractBetween(token,':O:',':R:');
R = extractBetween(token,':R:',':c:');
    s = strsplit(token,':c:');
c = s(2);
P = P{1};
C = C{1};
W = W{1};
F = F{1};
O = O{1};
R = R{1};
c = c{1};

%--------------------------------------------------------------------
function show_plot(hObject,handles)
%
    [s1, Q, param_ind] = select_data(handles,'Q1',2); % group of data - 1,2 for 1D histos, 3 is for 2D histos
    [s2, ~] = select_data(handles,'Q1',1);
    %       
    logscale = false;
    
    MEAN = [];
    STD = [];
    YLABEL = [];
    XLABEL = [];
    AXESMINMAX = [];
    
    AXES = handles.axes1;

    YLABEL = 'pdf';
    switch Q
        case 'roi.Ripley'
            YLABEL = 'AU';        
    end
    
    AXESMINMAX = [handles.minmaxlimits(param_ind,1) handles.minmaxlimits(param_ind,2)];
    
    switch Q
        case 'cl.MeanDoC'
            XLABEL = 'Degree of Co-localising';
        case 'cl.Area'
            logscale = true;
            AXESMINMAX = log10(AXESMINMAX);            
            XLABEL = 'log10(Area[{nm}^2])';
        case 'cl.Circularity'
            XLABEL = 'Circularity';
        case 'cl.Nb'
            logscale = true;
            AXESMINMAX = log10(AXESMINMAX);
            XLABEL = 'log10(N)';
        case 'cl.NormDensity'
            logscale = true;
            AXESMINMAX = log10(AXESMINMAX);
            XLABEL = 'log10(relative localisations density)';
        case 'cl.Elongation'
            XLABEL = 'Elongation';
        case 'cl.Density'
            logscale = true;
            AXESMINMAX = log10(AXESMINMAX);
            XLABEL = 'log10(localization density [1/{nm}^2)]';
        case 'roi.Ripley'
            XLABEL = 'distance [nm]';
            plot(AXES,handles.Ripley_distance,s1,'b.-',handles.Ripley_distance,s2,'r.-');
            if get(handles.legend,'Value')
                %legend(AXES,{[C1 ' ' Q ' channel ' num2str(channel1)],[C2 ' ' Q ' channel ' num2str(channel2)]});
            end
        case 'roi.Ripley(MAX)'            
            XLABEL = 'distance at RipleyK maximum [nm]';
        case 'roi.ClusterDensity'
            XLABEL = 'clusters density estimate [1/nm^2]';            
    end
                              
    if ~isempty(s1) && ~isempty(s2) && ~strcmp(Q,'roi.Ripley')
        
        s1 = s1(~isnan(s1));
        s2 = s2(~isnan(s2));
    
        MEAN1 = mean(s1);
        STD1= std(s1);
        MEAN2 = mean(s2);
        STD2 = std(s2);
                
        stat_vals = calculate_statistics(s1,s2,handles);
        set(handles.statistics_table,'Data',stat_vals);
        guidata(hObject, handles);
                
        if logscale
            s1 = log10(s1);
            s2 = log10(s2);
        end            
        
        % calculate stats - before log or after?
        %stat_vals = calculate_statistics(s1,s2);
        %set(handles.statistics_names,'Data',stat_vals); 
        %guidata(hObject, handles);
                        
        h1 = histogram(AXES,s1,'Normalization','pdf');
        hold(AXES,'on');
        h2 = histogram(AXES,s2,'Normalization','pdf');
        hold(AXES,'off');
            maxval = max(max(h1.Values),max(h2.Values));
                coef = 1;
        if ~isempty(AXESMINMAX)        
            axis(AXES,[AXESMINMAX 0. coef*maxval]);
        end
        
        if get(handles.legend,'Value')
            legend(AXES,{[Q ' mean: ' num2str(MEAN1) ' std: ' num2str(STD1)],[Q ' mean: ' num2str(MEAN2) ' std: ' num2str(STD2)]}, ...
            'FontSize',7.5);
        end 
    else
        set(handles.statistics_table,'Data',nan(5,1));
    end
        if get(handles.axes_log,'Value')
         set(AXES,'YScale','log'); % ?
        end
        
    xlabel(AXES,XLABEL,'fontsize',10);
    ylabel(AXES,YLABEL,'fontsize',10);
       
    grid(AXES,'on');
           
%--------------------------------------------------------------------
function [s, parameter,param_ind] = select_data(handles,table_token,group_index)
     table = eval(['handles.' table_token]); 
     %
     roi_data = handles.tot_data;
     %
     param_ind = get(table,'Value');
     parameter = char(handles.param_names(param_ind));
     %
     cluster_params = {'cl.Area','cl.Nb','cl.Circularity','cl.Elongation','cl.MeanDoC','cl.Density','cl.NormDensity'};     
     %
     if strcmp(parameter,'roi.Ripley')
         s = zeros(size(handles.Ripley_distance));
     elseif strcmp(parameter,'roi.Ripley(MAX)')
         s = nan(1,100000000,'single'); % enough to handle clusters...
     else %clusters
         s = nan(1,100000000,'single'); % enough to handle clusters...
     end
         tic
         cnt = 0;
             for k=1:numel(roi_data)
                    if ~handles.mask(k,group_index), continue, end
                    ass = handles.tot_data{k};
                    if isfield(ass,'DBSCAN_clusters') ...
                       && ~isempty(ass.DBSCAN_clusters) ...
                       && ismember(parameter,cluster_params)
                        %
                        if strcmp(parameter,'cl.NormDensity') % needs to calculate norma
                            [P,C,W,F,~,~,c] = get_attributes(ass.token);
                            index = handles.FOV_ind([':P:' P ':C:' C ':W:' W ':F:' F]);
                            cind = handles.Chan_ind(c);
                            norma = handles.FOVs_density(index,cind);
                        end
                        %
                            for cl = 1:numel(ass.DBSCAN_clusters)
                                cnt=cnt+1;
                                switch parameter
                                    case 'cl.Area'
                                        s(cnt) = ass.DBSCAN_clusters{cl}.Area;
                                    case 'cl.Nb'
                                        s(cnt) = ass.DBSCAN_clusters{cl}.Nb;
                                    case 'cl.Circularity'
                                        s(cnt) = ass.DBSCAN_clusters{cl}.Circularity;
                                    case 'cl.Elongation'
                                        s(cnt) = ass.DBSCAN_clusters{cl}.Elongation;
                                    case 'cl.MeanDoC'
                                        s(cnt) = ass.DBSCAN_clusters{cl}.MeanDoC;
                                    case 'cl.Density'
                                        s(cnt) = ass.DBSCAN_clusters{cl}.Nb/ass.DBSCAN_clusters{cl}.Area;
                                    case 'cl.NormDensity'
                                        s(cnt) = ass.DBSCAN_clusters{cl}.Nb/ass.DBSCAN_clusters{cl}.Area/norma;
                                end
                            end
                    end  
                    if strcmp(parameter,'roi.Ripley')
                        ripvals = ass.RipleyK_val;
                        if ~isempty(ripvals)
                            s = s + ripvals;
                            cnt = cnt + 1;
                        end
                    end
                    if strcmp(parameter,'roi.Ripley(MAX)')
                        maxval = handles.RipleyK_data_MAX(k);
                        if ~isnan(maxval)
                            cnt = cnt + 1;                            
                             s(cnt) = maxval;
                        end
                    end
                    if strcmp(parameter,'roi.ClusterDensity')
                        val = ass.ClusterDensity;
                        if ~isnan(val)
                            cnt = cnt + 1;                            
                             s(cnt) = val;
                        end
                    end
                    
             end        
         toc
     if ismember(parameter,cluster_params) || strcmp(parameter,'roi.Ripley(MAX)')
        s = s(1:cnt);
        numel(s);
     elseif strcmp(parameter,'roi.Ripley')
        s = s/cnt; 
     end

% %---------------------------------------------
function show_2d_histogram(handles)
% 
    % quantifiers
    s = get(handles.Q_X,'String');
    k = get(handles.Q_X,'Value');
    Q_X = s{k};
    s = get(handles.Q_Y,'String');
    k = get(handles.Q_Y,'Value');
    Q_Y = s{k};
    
    [sx, Q_X, param_ind] = select_data(handles,'Q_X',3); % group of data - 1,2 for 1D histos, 3 is for 2D histos
    [sy, Q_Y, ~] = select_data(handles,'Q_Y',3);
            
    YLABEL = [];
    XLABEL = [];        
            
    logscale_X = false;
    switch Q_X
        case 'cl.MeanDoC'
                XLABEL = 'Degree of Co-localising';
        case 'cl.Area'
                logscale_X = true;
                XLABEL = 'log10(Area[{nm}^2])';
        case 'cl.Circularity'
                XLABEL = 'Circularity';
        case 'cl.Elongation'
                XLABEL = 'Elongation';            
        case 'cl.Nb'
                logscale_X = true;
                XLABEL = 'log10(N)';
        case 'cl.Density'
                logscale_X = true;
                XLABEL = 'log10(localization density [1/{nm}^2])';
        case 'cl.NormDensity'
                logscale_X = true;
                XLABEL = 'log10(relative localisations density)';                
        case 'roi.Ripley(MAX)'
            XLABEL = 'distance at RipleyK maximum [nm]';
        case 'roi.ClusterDensity'            
            XLABEL = 'clusters density estimate [1/nm^2]';
            
    end
    
    logscale_Y = false;
    switch Q_Y
        case 'cl.MeanDoC'
                YLABEL = 'Degree of Co-localising';
        case 'cl.Area'
                logscale_Y = true;
                YLABEL = 'log10(Area[{nm}^2])';
        case 'cl.Circularity'
                YLABEL = 'Circularity';
        case 'cl.Elongation'
                YLABEL = 'Elongation';            
        case 'cl.Nb'
                logscale_Y = true;
                YLABEL = 'log10(N)';
        case 'cl.Density'
                logscale_Y = true;
                YLABEL = 'log10(localization density [1/{nm}^2])';
        case 'cl.NormDensity'
                logscale_Y = true;
                YLABEL = 'log10(relative localisations density)';                
        case 'roi.Ripley(MAX)'
            YLABEL = 'distance at RipleyK maximum [nm]';
        case 'roi.ClusterDensity'            
            YLABEL = 'clusters density estimate [1/nm^2]';            
    end
    %
    if ~isempty(sx) && ~isempty(sy) && length(sx)==length(sy)
        if logscale_X
            sx = log10(sx);
        end 
        if logscale_Y
            sy = log10(sy);
        end          

        nBins = 50;
        AXES = handles.corr_plot;
                             
        if ~get(handles.bar_histo,'Value')
            %histogram2(AXES,sx,sy,[nBins nBins],'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
            histogram2(AXES,sx,sy,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
        else
            histogram2(AXES,sx,sy,'DisplayStyle','bar3','ShowEmptyBins','on','Normalization','pdf','FaceColor','flat');                    
            if get(handles.bar_log,'Value')
             set(AXES,'zscale','log');
            end
        end
        
    xlabel(AXES,XLABEL,'fontsize',8);
    ylabel(AXES,YLABEL,'fontsize',8);         
        
    end    

%--------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
    
% % --- Executes on button press in axes_log.
function axes_log_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);

% % --- Executes on button press in remove_outliers.
function remove_outliers_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);

% % --- Executes on button press in legend.
function legend_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);

function stat_vals = calculate_statistics(x1,x2,handles)
    if get(handles.remove_outliers,'Value')
        x1=rmoutliers(x1,'median');
        x2=rmoutliers(x2,'median');                            
    end
    [~,P1] = kstest2(x1,x2);
    [~,P2] = ttest2(x1,x2,'vartype','unequal'); % not sure                        
    [P3,~] = ranksum(x1,x2);
    % Cohen's d
    N1 = numel(x1);
    N2 = numel(x2);
    s = sqrt( 1/(N1+N2)*( (N1-1)*var(x1) + (N2-1)*var(x2) ) );
    P4 = abs( mean(x1) - mean(x2) )/s;
    P5 = abs(median(x1)-median(x2));
stat_vals = [P1 P2 P3 P4 P5]';                        


% --------------------------------------------------------------------
function Tools_Callback(hObject, eventdata, handles)
% hObject    handle to Tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save_clusters_data_Callback(hObject, eventdata, handles)
     try
        roi_data = handles.tot_data;
     catch
        disp('no data');
     end
     if isempty(roi_data), return, end

    cluster_params = {'cl.Area','cl.Nb','cl.Circularity','cl.Elongation','cl.MeanDoC','cl.Density','cl.NormDensity'};
    is_Elongation = true;
    is_MeanDoC = true;
    if ~ismember(handles.param_names,'cl.Elongation')
        setxor(cluster_params,'cl.Elongation');
        is_Elongation = false;
    end
    if ~ismember(handles.param_names,'cl.MeanDoC')
        setxor(cluster_params,'cl.MeanDoC');
        is_MeanDoC = false;
    end        
                            %
                            captions = {'Plate','Condition','Well','FOV','Object','ROI','channel','Xc','Yc'};
                            captions = [captions 'Area'];
                            captions = [captions 'Nb'];
                            captions = [captions 'Circularity'];
                                        if is_Elongation
                                            captions = [captions 'Elongation'];
                                        end
                                        if is_MeanDoC
                                            captions = [captions 'MeanDoC'];
                                        end
                            captions = [captions 'loc.Density'];
                            captions = [captions 'norm loc.Density'];   
    
         row_size = numel(captions);
             
         s = cell(100000000,row_size); % enough to handle clusters...
         cnt = 0;
         tic
         w = waitbar(0, 'grabbing clusters...');
         N_rois = numel(roi_data);
             for k=1:N_rois 
                 waitbar(k/N_rois,w);
                    ass = roi_data{k};
                    if isfield(ass,'DBSCAN_clusters') && ~isempty(ass.DBSCAN_clusters)
                        %
                            [P,C,W,F,O,R,c] = get_attributes(ass.token);
                            index = handles.FOV_ind([':P:' P ':C:' C ':W:' W ':F:' F]);
                            cind = handles.Chan_ind(c);
                            norma = handles.FOVs_density(index,cind);
                            
                            for cl = 1:numel(ass.DBSCAN_clusters)
                                cnt=cnt+1;
                                        s{cnt,1} = P;
                                        s{cnt,2} = C;
                                        s{cnt,3} = W;
                                        s{cnt,4} = F;
                                        s{cnt,5} = O;
                                        s{cnt,6} = R;
                                        s{cnt,7} = c;                                        
                                        p=1;
                                        s{cnt,7+p} = ass.DBSCAN_clusters{cl}.Xc;
                                        p=p+1;
                                        s{cnt,7+p} = ass.DBSCAN_clusters{cl}.Yc;                                                                                                                        
                                        p=p+1;
                                        s{cnt,7+p} = ass.DBSCAN_clusters{cl}.Area;
                                        p=p+1;
                                        s{cnt,7+p} = ass.DBSCAN_clusters{cl}.Nb;
                                        p=p+1;                                        
                                        s{cnt,7+p} = ass.DBSCAN_clusters{cl}.Circularity;
                                        p=p+1;
                                        if is_Elongation
                                            s{cnt,7+p} = ass.DBSCAN_clusters{cl}.Elongation;                                        
                                            p=p+1;                                        
                                        end
                                        if is_MeanDoC
                                            s{cnt,7+p} = ass.DBSCAN_clusters{cl}.MeanDoC;                                                   
                                            p=p+1;                                        
                                        end
                                        s{cnt,7+p} = ass.DBSCAN_clusters{cl}.Nb/ass.DBSCAN_clusters{cl}.Area;
                                        p=p+1;                                        
                                        s{cnt,7+p} = ass.DBSCAN_clusters{cl}.Nb/ass.DBSCAN_clusters{cl}.Area/norma;
                            end
                    end  
             end  
             waitbar(1, w);
             close(w);             
         s = s(1:cnt,:);
                                             
        [filename, pathname] = uiputfile('*.csv', 'save clusters data as');
        fullfilename = [pathname filesep filename];
        s = [captions; s];
        cell2csv(fullfilename,s); % very slow
        toc
                
        
function ret = get_cluster_density_per_ROI(handles,Lmax)
ret = cell(size(handles.tot_data));
for k=1:numel(handles.tot_data)
    %
    ass = handles.tot_data{k};
    ass.ClusterDensity = nan;
    ret{k} = ass;    
    try % to calculate clusters density
    if isfield(ass,'DBSCAN_clusters') && ~isempty(ass.DBSCAN_clusters)
        % 
        points = zeros(numel(ass.DBSCAN_clusters),2); % set of points
        for cl = 1:numel(ass.DBSCAN_clusters)
            points(cl,:) = [ass.DBSCAN_clusters{cl}.Xc ass.DBSCAN_clusters{cl}.Yc];
        end

         dt = delaunayTriangulation(points); 
         if isempty(dt), continue, end
%         tri = dt.ConnectivityList;
%         g = digraph(tri, tri(:, [2 3 1]));
%         A = adjacency(g);
%         A = A | A';
%         g = graph(A);
%         figure(22);        
%             tri = dt.ConnectivityList;
%             g = digraph(tri, tri(:, [2 3 1]));
%             A = adjacency(g);
%             A = A | A';
%             g = graph(A);
%             plot(g);        
%          figure(23);
%              IC = incenter(dt);
%              triplot(dt)
%              hold(gca,'on');
%              plot(IC(:,1),IC(:,2),'*r')        
%              hold(gca,'off');                                                                            
         
         EDGES = edges(dt);
%         figure(22)
%         cla(gca,'reset');
         acc = [];
         for ei=1:size(EDGES,1)
             v1 = EDGES(ei,1);
             v2 = EDGES(ei,2);
             p1 = dt.Points(v1,:);
             p2 = dt.Points(v2,:);
             li = norm(p1-p2);             
             if li<Lmax
%                 hold('on');
%             line(gca,[p1(1) p2(1)],[p1(2) p2(2)],'color','b');
%                 hold('off');
                 acc = [acc; li];
             else
%                  hold('on');
%              line([p1(1) p2(1)],[p1(2) p2(2)],'color','r');
%                  hold('off');
             end
         end
         if ~isempty(acc)
            ass.ClusterDensity = 2/mean(acc)^2;
            ret{k} = ass;
         end         
    end
    catch
    end
end








% --- Executes on selection change in Q2_Y.
function Q2_Y_Callback(hObject, eventdata, handles)
show_2d_histogram2(handles); 


% --- Executes during object creation, after setting all properties.
function Q2_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q2_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Q2_X.
function Q2_X_Callback(hObject, eventdata, handles)
show_2d_histogram2(handles); 

% --- Executes during object creation, after setting all properties.
function Q2_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q2_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bar2_histo.
function bar2_histo_Callback(hObject, eventdata, handles)
show_2d_histogram2(handles); 

% --- Executes on button press in bar2_log.
function bar2_log_Callback(hObject, eventdata, handles)
show_2d_histogram2(handles); 

% Hint: get(hObject,'Value') returns toggle state of bar2_log


% %---------------------------------------------
function show_2d_histogram2(handles)
% 
    % quantifiers
    s = get(handles.Q2_X,'String');
    k = get(handles.Q2_X,'Value');
    Q_X = s{k};
    s = get(handles.Q2_Y,'String');
    k = get(handles.Q2_Y,'Value');
    Q_Y = s{k};
    
    [sx, Q_X, param_ind] = select_data(handles,'Q2_X',4); % group of data - 1,2 for 1D histos, 3 is for 2D histos
    [sy, Q_Y, ~] = select_data(handles,'Q2_Y',4);
            
    YLABEL = [];
    XLABEL = [];        
            
    logscale_X = false;
    switch Q_X
        case 'cl.MeanDoC'
                XLABEL = 'Degree of Co-localising';
        case 'cl.Area'
                logscale_X = true;
                XLABEL = 'log10(Area[{nm}^2])';
        case 'cl.Circularity'
                XLABEL = 'Circularity';
        case 'cl.Elongation'
                XLABEL = 'Elongation';            
        case 'cl.Nb'
                logscale_X = true;
                XLABEL = 'log10(N)';
        case 'cl.Density'
                logscale_X = true;
                XLABEL = 'log10(localization density [1/{nm}^2])';
        case 'cl.NormDensity'
                logscale_X = true;
                XLABEL = 'log10(relative localisations density)';                
        case 'roi.Ripley(MAX)'
            XLABEL = 'distance at RipleyK maximum [nm]';
        case 'roi.ClusterDensity'            
            XLABEL = 'clusters density estimate [1/nm^2]';
            
    end
    
    logscale_Y = false;
    switch Q_Y
        case 'cl.MeanDoC'
                YLABEL = 'Degree of Co-localising';
        case 'cl.Area'
                logscale_Y = true;
                YLABEL = 'log10(Area[{nm}^2])';
        case 'cl.Circularity'
                YLABEL = 'Circularity';
        case 'cl.Elongation'
                YLABEL = 'Elongation';            
        case 'cl.Nb'
                logscale_Y = true;
                YLABEL = 'log10(N)';
        case 'cl.Density'
                logscale_Y = true;
                YLABEL = 'log10(localization density [1/{nm}^2])';
        case 'cl.NormDensity'
                logscale_Y = true;
                YLABEL = 'log10(relative localisations density)';                
        case 'roi.Ripley(MAX)'
            YLABEL = 'distance at RipleyK maximum [nm]';
        case 'roi.ClusterDensity'            
            YLABEL = 'clusters density estimate [1/nm^2]';            
    end
    %
    if ~isempty(sx) && ~isempty(sy) && length(sx)==length(sy)
        if logscale_X
            sx = log10(sx);
        end 
        if logscale_Y
            sy = log10(sy);
        end          

        nBins = 50;
        AXES = handles.corr_plot2;
                             
        if ~get(handles.bar2_histo,'Value')
            %histogram2(AXES,sx,sy,[nBins nBins],'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
            histogram2(AXES,sx,sy,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
        else
            histogram2(AXES,sx,sy,'DisplayStyle','bar3','ShowEmptyBins','on','Normalization','pdf','FaceColor','flat');                    
            if get(handles.bar_log,'Value')
             set(AXES,'zscale','log');
            end
        end
        
    xlabel(AXES,XLABEL,'fontsize',8);
    ylabel(AXES,YLABEL,'fontsize',8);         
        
    end    