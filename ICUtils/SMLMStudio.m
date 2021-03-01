function varargout = SMLMStudio(varargin)
% SMLMSTUDIO MATLAB code for SMLMStudio.fig
%      SMLMSTUDIO, by itself, creates a new SMLMSTUDIO or raises the existing
%      singleton*.
%
%      H = SMLMSTUDIO returns the handle to a new SMLMSTUDIO or the handle to
%      the existing singleton*.
%
%      SMLMSTUDIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMLMSTUDIO.M with the given input arguments.
%
%      SMLMSTUDIO('Property','Value',...) creates a new SMLMSTUDIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SMLMStudio_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SMLMStudio_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SMLMStudio

% Last Modified by GUIDE v2.5 01-Mar-2021 15:58:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SMLMStudio_OpeningFcn, ...
                   'gui_OutputFcn',  @SMLMStudio_OutputFcn, ...
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


% --- Executes just before SMLMStudio is made visible.
function SMLMStudio_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SMLMStudio (see VARARGIN)

% Choose default command line output for SMLMStudio
handles.output = hObject;

handles.figureName = get(handles.figure1,'Name');

handles.data = [];
handles.filenames = [];
%set(handles.chart_statistic,'String',{'p KS','p t-test','p Wilcoxon','Cohen d'});

set(handles.Ch1,'String',{'462','647'});
set(handles.Ch2,'String',{'462','647'});
set(handles.Ch_Y,'String',{'462','647'});

handles.ROI_side = 2000; % [nm], square ROI

handles.param_names = {'cl.DoC','cl.Area','cl.Circularity','cl.Nloc','cl.Density','roi.Ripley','roi.RDF','roi.xRipley','roi.Ripley(MAX)','roi.RDF(MAX)','roi.xRipley(MAX)','roi.SAA','roi.SAA(ratio)'};

set(handles.Q1,'String',handles.param_names);
set(handles.C1,'String',{'H3K27me3','H3K9ac','H3K9me3'});
set(handles.C2,'String',{'H3K27me3','H3K9ac','H3K9me3'});
set(handles.Q_X,'String',handles.param_names);
set(handles.Q_Y,'String',handles.param_names);
set(handles.C_Y,'String',{'H3K27me3','H3K9ac','H3K9me3'});

set(handles.minmax_table, 'RowName',handles.param_names);
set(handles.minmax_table, 'ColumnName', {'min','max'});
%set(handles.minmax_table,'CellEditCallback',@minmax_check_callback);

handles.statistics_names = {'p-value: KS','p-value: t-test','p-value: Wilcoxon','Cohen"s d','|median diff|'};
set(handles.statistics_table, 'RowName',handles.statistics_names);
set(handles.statistics_table, 'ColumnName', {'value'});

set(handles.bar_log,'Enable','off');

handles.DOC_DBSCAN_samples = [];

handles.minmaxlimits = define_minmax(handles);
set(handles.minmax_table,'Data',handles.minmaxlimits);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SMLMStudio wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SMLMStudio_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in C1.
function C1_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);

% --- Executes during object creation, after setting all properties.
function C1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Q1.
function Q1_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Q1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ch1.
function Ch1_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Ch1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in C2.
function C2_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);


% --- Executes during object creation, after setting all properties.
function C2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Q2.
function Q2_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Q2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ch2.
function Ch2_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Ch2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% % --------------------------------------------------------------------
function load_data_Callback(hObject, eventdata, handles)
[filenames,pathname] = uigetfile('*.mat','Select track data files',pwd,'MultiSelect','on');                
if isempty(filenames), return, end       
if isnumeric(filenames) && 0==filenames, return, end

handles.filenames = filenames;
handles.conditions_Ch1 = cell(0);
handles.conditions_Ch2 = cell(0);
handles.conditions = cell(0);

handles.data = cell(0);
for k=1:numel(handles.filenames)
    load([pathname filesep handles.filenames{k}]);
    handles.data{k} = data;
    handles.conditions_Ch1{k} = data.Ch1_label;    
    handles.conditions_Ch2{k} = data.Ch2_label;
    k
end

UC1 = unique(handles.conditions_Ch1);
UC2 = unique(handles.conditions_Ch2);
if numel(UC1)>1
    UC = UC1;
    handles.conditions = handles.conditions_Ch1;
elseif numel(UC2)>1
    UC = UC2;    
    handles.conditions = handles.conditions_Ch2;    
end
set(handles.C1,'String',UC);
set(handles.C2,'String',UC);
    
handles.DOC_DBSCAN_samples = get_DBSCAN_DoC_samples(handles);

[handles.RipleyK_x, ...
 handles.RipleyK_curves, ...
 handles.RipleyK_curves_std, ...
 handles.RipleyK_max_sample] = get_Ripley_samples(handles);

[handles.RDF_x, ...
 handles.RDF_curves, ...
 handles.RDF_curves_std, ...
 handles.RDF_max_sample] = get_RDF_samples(handles);

[handles.Hr_x, ...
 handles.Hr_curves, ...
 handles.Hr_curves_std, ...
 handles.Hr_max_sample] = get_Hr_samples(handles);

[handles.SAA_x, ...
 handles.SAA_curves, ...
 handles.SAA_curves_std, ...
 handles.sim_SAA_curves, ...
 handles.sim_SAA_curves_std, ...
 handles.SAA_fraction_sample, ...
 handles.SAA_cutoff] = get_SAA_samples(handles);

handles.minmaxlimits = define_minmax(handles);
set(handles.minmax_table,'Data',handles.minmaxlimits);

% Update handles structure
guidata(hObject, handles);

set(handles.figure1, 'Name', [handles.figureName ' : ' pathname]);

show_plot(hObject,handles);
show_2d_histogram(handles);

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function sample = get_DBSCAN_DoC_samples(handles)
%
%1 Area
%2 Circularity
%3 DoC
%4 Nloc
%5 TO DO - relative density 

UC = get(handles.C1,'String');

sample = cell(numel(UC),2); % 2 channels

for k=1:numel(handles.data)
        index = find(ismember(UC,handles.conditions{k}));

        sample_cur1 = handles.data{k}.DoC_DBSCAN_clusters_ch1; % first channel
        sample_cur2 = handles.data{k}.DoC_DBSCAN_clusters_ch2; % second channel
        
        sample_cur1(sample_cur1(:,2)>1,2)=1; % fix Circularity
        sample_cur2(sample_cur2(:,2)>1,2)=1;
                
        if 4 == size(sample_cur1,2) % SHOULD BE 5 - made on compacting
%             % density involution % INAPROPRIATE!!!
%             nROIs = max(length(handles.data{k}.SAAstruct),length(handles.data{k}.SpatialData_Ch1));
%             sROIS = nROIs*(handles.ROI_side)^2; % nm
%             density_k_1 = sum(sample_cur1(:,4))/sROIS;
%             rel_density_1 = sample_cur1(:,4)./sample_cur1(:,1)/density_k_1;
%             density_k_2 = sum(sample_cur2(:,4))/sROIS;
%             rel_density_2 = sample_cur2(:,4)./sample_cur2(:,1)/density_k_2;

            density_k_1 = sum(sample_cur1(:,4))/sum(sample_cur1(:,1));
            rel_density_1 = sample_cur1(:,4)./sample_cur1(:,1)/density_k_1;
            density_k_2 = sum(sample_cur2(:,4))/sum(sample_cur2(:,1));
            rel_density_2 = sample_cur2(:,4)./sample_cur2(:,1)/density_k_2;
        end
                
        sample{index,1} = [sample{index,1}; [sample_cur1 rel_density_1]]; % append
        sample{index,2} = [sample{index,2}; [sample_cur2 rel_density_2]];       
end

%--------------------------------------------------------
function [xaxis,curves,curves_std,max_sample] = get_Ripley_samples(handles)
%
UC = get(handles.C1,'String');

max_sample = cell(numel(UC),2); % 2 channels
curves = cell(numel(UC),2);
curves_std = cell(numel(UC),2);

xaxis = handles.data{1}.RipelyKx;

curves_acc = cell(numel(UC),2);

% accumulate
for k=1:numel(handles.data)
        cond_index = find(ismember(UC,handles.conditions{k}));
            
        curves1 = handles.data{k}.RipelyK_Ch1_data;
        curves2 = handles.data{k}.RipelyK_Ch2_data;
        
        % normalize - max, mean, or median?
        %curves1 = curves1./max(curves1,[],1);
        %curves2 = curves2./max(curves2,[],1);
        
        curves_acc{cond_index,1} = [curves_acc{cond_index,1} curves1];
        curves_acc{cond_index,2} = [curves_acc{cond_index,2} curves2];        
end
%
for k=1:numel(UC) % over conditions
        curves{k,1} = median(curves_acc{k,1},2);
        curves{k,2} = median(curves_acc{k,2},2);
        curves_std{k,1} = std(curves_acc{k,1},0,2);
        curves_std{k,2} = std(curves_acc{k,2},0,2);    
end

for k=1:numel(UC) % over conditions
    for c=1:2
        curves_kc = curves_acc{k,c};
        for r=1:size(curves_kc,2) % over ROIs
            curv_r = curves_kc(:,r);
            maxind = find(curv_r==max(curv_r));
            maxind = maxind(1); % :)
            max_sample{k,c} = [max_sample{k,c}; xaxis(maxind)];
        end
    end
end


%--------------------------------------------------------
function show_plot(hObject,handles) % quantifier, condition

    % quantifier
    s = get(handles.Q1,'String');
    k = get(handles.Q1,'Value');
    Q = s{k};    
    % channel
    channel1 = get(handles.Ch1,'Value');
    channel2 = get(handles.Ch2,'Value');
    % condition
         k = get(handles.C1,'Value');
         s = get(handles.C1,'String');
         C1 = s{k};    
         k = get(handles.C2,'Value');
         s = get(handles.C2,'String');
         C2 = s{k};        
    UC = get(handles.C1,'String');
    index1 = find(ismember(UC,C1));
    index2 = find(ismember(UC,C2));

    s1 = []; % samples for histograms
    s2 = [];
        
    logscale = false;
    
    MEAN = [];
    STD = [];
    YLABEL = [];
    XLABEL = [];
    AXESMINMAX = [];
    
    AXES = handles.axes1;

    switch Q
        case {'cl.DoC','cl.Area','cl.Circularity','cl.Nloc'}
            YLABEL = 'pdf';        
    end
    
    switch Q
        case 'cl.DoC'
            s1 = handles.DOC_DBSCAN_samples{index1,channel1}(:,3);
            s2 = handles.DOC_DBSCAN_samples{index2,channel2}(:,3);
            XLABEL = 'Degree of Co-localising';
            AXESMINMAX = [handles.minmaxlimits(1,1) handles.minmaxlimits(1,2)];
        case 'cl.Area'
            s1 = handles.DOC_DBSCAN_samples{index1,channel1}(:,1);
            s2 = handles.DOC_DBSCAN_samples{index2,channel2}(:,1);
            XLABEL = 'log10(Area[nm])';
            AXESMINMAX = log10([handles.minmaxlimits(2,1) handles.minmaxlimits(2,2)]);
                logscale = true;
        case 'cl.Circularity'
            s1 = handles.DOC_DBSCAN_samples{index1,channel1}(:,2);
            s2 = handles.DOC_DBSCAN_samples{index2,channel2}(:,2);
            XLABEL = 'Circularity';
            AXESMINMAX = [handles.minmaxlimits(3,1) handles.minmaxlimits(3,2)];
        case 'cl.Nloc'
            s1 = handles.DOC_DBSCAN_samples{index1,channel1}(:,4);
            s2 = handles.DOC_DBSCAN_samples{index2,channel2}(:,4);
                logscale = true;
            XLABEL = 'log10(N)';
            AXESMINMAX = log10([handles.minmaxlimits(4,1) handles.minmaxlimits(4,2)]);
        case 'cl.Density'
            s1 = handles.DOC_DBSCAN_samples{index1,channel1}(:,5);
            s2 = handles.DOC_DBSCAN_samples{index2,channel2}(:,5);
                logscale = true;
            XLABEL = 'log10(relative localisations density)';
            AXESMINMAX = log10([handles.minmaxlimits(5,1) handles.minmaxlimits(5,2)]);            
        case 'roi.Ripley'
            c1 = handles.RipleyK_curves{index1,channel1};
            c2 = handles.RipleyK_curves{index2,channel2};
            XLABEL = 'distance [nm]';
            YLABEL = 'AU';
            plot(AXES,handles.RipleyK_x,c1,'b.-',handles.RipleyK_x,c2,'r.-');
            if get(handles.legend,'Value')
                legend(AXES,{[C1 ' ' Q ' channel ' num2str(channel1)],[C2 ' ' Q ' channel ' num2str(channel2)]});
            end
        case 'roi.RDF'
            c1 = handles.RDF_curves{index1,channel1};
            c2 = handles.RDF_curves{index2,channel2};
            XLABEL = 'distance [nm]';
            YLABEL = 'AU';
            plot(AXES,handles.RDF_x,c1,'b.-',handles.RDF_x,c2,'r.-');
            if get(handles.legend,'Value')
                legend(AXES,{[C1 ' ' Q ' channel ' num2str(channel1)],[C2 ' ' Q ' channel ' num2str(channel2)]});
            end
        case 'roi.xRipley'
            c1 = handles.Hr_curves{index1,channel1};
            c2 = handles.Hr_curves{index2,channel2};
            XLABEL = 'distance [nm]';
            YLABEL = 'AU';
            plot(AXES,handles.Hr_x,c1,'b.-',handles.Hr_x,c2,'r.-');
            if get(handles.legend,'Value')
                legend(AXES,{[C1 ' ' Q ' channel ' num2str(channel1)],[C2 ' ' Q ' channel ' num2str(channel2)]});
            end            
        case 'roi.Ripley(MAX)'            
            s1 = handles.RipleyK_max_sample{index1,channel1};
            s2 = handles.RipleyK_max_sample{index2,channel2};
            XLABEL = 'distance at RipleyK maximum [nm]';
            AXESMINMAX = [handles.minmaxlimits(9,1) handles.minmaxlimits(9,2)];
        case 'roi.RDF(MAX)'
            s1 = handles.RDF_max_sample{index1,channel1};
            s2 = handles.RDF_max_sample{index2,channel2};
            XLABEL = 'distance at RDF maximum [nm]';
            AXESMINMAX = []; %[handles.minmaxlimits(9,1) handles.minmaxlimits(9,2)];            
        case 'roi.xRipley(MAX)'
            s1 = handles.Hr_max_sample{index1,channel1};
            s2 = handles.Hr_max_sample{index2,channel2};
            XLABEL = 'distance at cross-Ripley maximum [nm]';
            AXESMINMAX = []; %[handles.minmaxlimits(9,1) handles.minmaxlimits(9,2)];                        
        case 'roi.SAA'
            c1 = handles.SAA_curves{index1,channel1};
            c2 = handles.SAA_curves{index2,channel2};
            sim_c1 = handles.sim_SAA_curves{index1,channel1};
            sim_c2 = handles.sim_SAA_curves{index2,channel2};
            %            
            XLABEL = 'distance [nm]';
            YLABEL = 'AU';
            x_cutoff1 = mean(handles.SAA_cutoff{index1});
            x_cutoff2 = mean(handles.SAA_cutoff{index2});
            s = [c1 c2 sim_c1 sim_c2];
            ymax = max(s(:));
            plot(AXES,handles.SAA_x,c1,'b.-',handles.SAA_x,c2,'r.-',handles.SAA_x,sim_c1,'b:',handles.SAA_x,sim_c2,'r:','linewidth',1.5);
            hold(AXES,'on');
            plot(AXES,[x_cutoff1 x_cutoff1],[1e-10 ymax],'b--',[x_cutoff2 x_cutoff2],[1e-10 ymax],'r--');
            hold(AXES,'off');            
            if get(handles.legend,'Value')
                                
                s1_ = handles.SAA_fraction_sample{index1,channel1};
                s2_ = handles.SAA_fraction_sample{index2,channel2};
                f1 = mean(s1_);
                f2 = mean(s2_);
                n1 = numel(s1_);
                n2 = numel(s2_);                
                legend(AXES,{[C1 ' ' Q ' channel ' num2str(channel1) ', N= ' num2str(n1) ', f_{avr}= ' num2str(f1)], ... 
                    [C2 ' ' Q ' channel ' num2str(channel2) ', N= ' num2str(n2) ', f_{avr}= ' num2str(f2)], ...
                    'sim','sim','cutoff','cutoff'},'FontSize',7.5);
                
            end            
        case 'roi.SAA(ratio)'
            s1 = handles.SAA_fraction_sample{index1,channel1};
            s2 = handles.SAA_fraction_sample{index2,channel2};
            XLABEL = 'interacting molecules fraction';
            AXESMINMAX = []; %[handles.minmaxlimits(9,1) handles.minmaxlimits(9,2)];
    end
        
    MEAN1 = mean(s1);
    STD1= std(s1);
    MEAN2 = mean(s2);
    STD2 = std(s2);
            
    if ~isempty(s1) && ~isempty(s2)
        
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
            legend(AXES,{[C1 ' ' Q ' mean: ' num2str(MEAN1) ' std: ' num2str(STD1)],[C2 ' ' Q ' mean: ' num2str(MEAN2) ' std: ' num2str(STD2)]}, ...
            'FontSize',7.5);
        end        
    end
        if get(handles.axes_log,'Value')
         set(AXES,'YScale','log'); % ?
        end
        
    xlabel(AXES,XLABEL,'fontsize',10);
    ylabel(AXES,YLABEL,'fontsize',10);
       
    grid(AXES,'on');
            
             
% --- Executes on button press in axes_log.
function axes_log_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);

% --- Executes on button press in remove_outliers.
function remove_outliers_Callback(hObject, eventdata, handles)
show_plot(hObject,handles);


% --- Executes on selection change in Q_Y.
function Q_Y_Callback(hObject, eventdata, handles)
show_2d_histogram(handles); 


% --- Executes during object creation, after setting all properties.
function Q_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in C_Y.
function C_Y_Callback(hObject, eventdata, handles)
show_2d_histogram(handles); 

% --- Executes during object creation, after setting all properties.
function C_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ch_Y.
function Ch_Y_Callback(hObject, eventdata, handles)
show_2d_histogram(handles); 


% --- Executes during object creation, after setting all properties.
function Ch_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ch_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Q_X.
function Q_X_Callback(hObject, eventdata, handles)
show_2d_histogram(handles); 

% --- Executes during object creation, after setting all properties.
function Q_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%---------------------------------------------
function show_2d_histogram(handles)

    % quantifiers
    s = get(handles.Q_X,'String');
    k = get(handles.Q_X,'Value');
    Q_X = s{k};
    s = get(handles.Q_Y,'String');
    k = get(handles.Q_Y,'Value');
    Q_Y = s{k};
    
    % channel
    channel_Y = get(handles.Ch_Y,'Value');
    % condition
         k = get(handles.C_Y,'Value');
         s = get(handles.C_Y,'String');
         C_Y = s{k};        
    UC = get(handles.C1,'String');
    index_Y = find(ismember(UC,C_Y));
    
    YLABEL = [];
    XLABEL = [];        
    
    sx = [];
    sy = [];
    
    logscale_X = false;
    switch Q_X
        case 'cl.DoC'
            sx = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,3);
            XLABEL = 'Degree of Co-localising';
        case 'cl.Area'
            sx = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,1);
                logscale_X = true;
                XLABEL = 'log10(Area[nm])';
        case 'cl.Circularity'
            sx = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,2);
            XLABEL = 'Circularity';
        case 'cl.Nloc'
            sx = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,4);            
                logscale_X = true;
                XLABEL = 'log10(N)';
        case 'cl.Density'
            sx = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,5);            
                logscale_X = true;
                XLABEL = 'log10(relative localisations density)';
        case 'roi.Ripley'
        case 'roi.RDF'
        case 'roi.xRipley'
        case 'roi.Ripley(MAX)'
            sx = handles.RipleyK_max_sample{index_Y,channel_Y};
            XLABEL = 'distance at RipleyK maximum [nm]';            
        case 'roi.RDF(MAX)'
        case 'roi.xRipley(MAX)'
        case 'roi.SAA'
        case 'roi.SAA(ratio)'
            sx = handles.SAA_fraction_sample{index_Y,channel_Y};            
            XLABEL = 'interacting molecules fraction';            
    end
    
    logscale_Y = false;
    switch Q_Y
        case 'cl.DoC'
            sy = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,3);
            YLABEL = 'Degree of Co-localising';
        case 'cl.Area'
            sy = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,1);
                logscale_Y = true;
                YLABEL = 'log10(Area[nm])';
        case 'cl.Circularity'
            sy = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,2);
            YLABEL = 'Circularity';
        case 'cl.Nloc'
            sy = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,4);            
                logscale_Y = true;
                YLABEL = 'log10(N)';
        case 'cl.Density'
            sy = handles.DOC_DBSCAN_samples{index_Y,channel_Y}(:,5);            
                logscale_Y = true;
                YLABEL = 'log10(relative localisations density)';
        case 'roi.Ripley'
        case 'roi.RDF'
        case 'roi.xRipley'
        case 'roi.Ripley(MAX)'
            sy = handles.RipleyK_max_sample{index_Y,channel_Y};
            YLABEL = 'distance at RipleyK maximum [nm]';                        
        case 'roi.RDF(MAX)'
        case 'roi.xRipley(MAX)'
        case 'roi.SAA'
        case 'roi.SAA(ratio)'
            sy = handles.SAA_fraction_sample{index_Y,channel_Y};            
            YLABEL = 'interacting molecules fraction';
    end
    
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
        

    
function minmaxlimits = define_minmax(handles)
minmaxlimits = zeros(numel(handles.param_names),2);

    if ~isempty(handles.DOC_DBSCAN_samples)

        s = get(handles.C2,'String');
        max_cond_ind = length(s);
        DoC = [];
        Area = [];
        Circularity = [];
        N = [];
        relDensity = [];

        for k=1:max_cond_ind
            for m=1:2
                DoC = [DoC; handles.DOC_DBSCAN_samples{k,m}(:,3)];
                Area = [Area; handles.DOC_DBSCAN_samples{k,m}(:,1)];
                Circularity = [Circularity; handles.DOC_DBSCAN_samples{k,m}(:,2)];
                N = [N; handles.DOC_DBSCAN_samples{k,m}(:,4)];
                relDensity = [relDensity; handles.DOC_DBSCAN_samples{k,m}(:,5)];                
            end
        end
        minmaxlimits(1,1)=-1;%min(DoC);
        minmaxlimits(1,2)=1;%max(DoC);
            minmaxlimits(2,1)=min(Area);
            minmaxlimits(2,2)=max(Area);
                minmaxlimits(3,1)=0; %min(Circularity);
                minmaxlimits(3,2)=1; %max(Circularity);
                    minmaxlimits(4,1)=3; %min(N);
                    minmaxlimits(4,2)=1e5; %max(N);
                        minmaxlimits(5,1)=1e-2; %min(relDensity);
                        minmaxlimits(5,2)=1e3; %max(relDensity);        
                            minmaxlimits(9,1)=0; %min(RipleyMaxdist);
                            minmaxlimits(9,2)=1000; %max(RipleyMaxdist);        
                        
    else
        return;
    end


% --- Executes on button press in legend.
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


%--------------------------------------------------------
function [xaxis,curves,curves_std,max_sample] = get_RDF_samples(handles)

UC = get(handles.C1,'String');

max_sample = cell(numel(UC),2); % 2 channels
curves = cell(numel(UC),2);
curves_std = cell(numel(UC),2);

curves_acc = cell(numel(UC),2);

xaxis = handles.data{1}.SpatialData_Ch1{1}.Xr;


% accumulate
for k=1:numel(handles.data)
        cond_index = find(ismember(UC,handles.conditions{k}));
 
        for m=1:numel(handles.data{k}.SpatialData_Ch1)
            if ~isempty(handles.data{k}.SpatialData_Ch1{m}) && ~isempty(handles.data{k}.SpatialData_Ch2{m})
                curve1 = handles.data{k}.SpatialData_Ch1{m}.Gr;
                curve2 = handles.data{k}.SpatialData_Ch2{m}.Gr;

                % normalize - max, mean, or median?
                %curves1 = curves1./max(curves1,[],1);
                %curves2 = curves2./max(curves2,[],1);

                curves_acc{cond_index,1} = [curves_acc{cond_index,1} curve1];
                curves_acc{cond_index,2} = [curves_acc{cond_index,2} curve2]; 
            end
        end
end
%
for k=1:numel(UC) % over conditions
        curves{k,1} = median(curves_acc{k,1},2);
        curves{k,2} = median(curves_acc{k,2},2);
        curves_std{k,1} = std(curves_acc{k,1},0,2);
        curves_std{k,2} = std(curves_acc{k,2},0,2);    
end

for k=1:numel(UC) % over conditions
    for c=1:2
        curves_kc = curves_acc{k,c};
        for r=1:size(curves_kc,2) % over ROIs
            curv_r = curves_kc(:,r);
            maxind = find(curv_r==max(curv_r));
            maxind = maxind(1); % :)
            max_sample{k,c} = [max_sample{k,c}; xaxis(maxind)];
        end
    end
end

%------------ this is full copy of the function above..
function [xaxis,curves,curves_std,max_sample] = get_Hr_samples(handles)

UC = get(handles.C1,'String');

max_sample = cell(numel(UC),2); % 2 channels
curves = cell(numel(UC),2);
curves_std = cell(numel(UC),2);

curves_acc = cell(numel(UC),2);

xaxis = handles.data{1}.SpatialData_Ch1{1}.Xr;

% accumulate
for k=1:numel(handles.data)
        cond_index = find(ismember(UC,handles.conditions{k}));
 
        for m=1:numel(handles.data{k}.SpatialData_Ch1)
            if ~isempty(handles.data{k}.SpatialData_Ch1{m}) && ~isempty(handles.data{k}.SpatialData_Ch2{m})
                curve1 = handles.data{k}.SpatialData_Ch1{m}.Hr(:,1); % ????
                curve2 = handles.data{k}.SpatialData_Ch2{m}.Hr(:,1);

                % normalize - max, mean, or median?
                %curves1 = curves1./max(curves1,[],1);
                %curves2 = curves2./max(curves2,[],1);

                curves_acc{cond_index,1} = [curves_acc{cond_index,1} curve1];
                curves_acc{cond_index,2} = [curves_acc{cond_index,2} curve2]; 
            end
        end
end
%
for k=1:numel(UC) % over conditions
        curves{k,1} = median(curves_acc{k,1},2);
        curves{k,2} = median(curves_acc{k,2},2);
        curves_std{k,1} = std(curves_acc{k,1},0,2);
        curves_std{k,2} = std(curves_acc{k,2},0,2);    
end

for k=1:numel(UC) % over conditions
    for c=1:2
        curves_kc = curves_acc{k,c};
        for r=1:size(curves_kc,2) % over ROIs
            curv_r = curves_kc(:,r);
            maxind = find(curv_r==max(curv_r));
            maxind = maxind(1); % :)
            max_sample{k,c} = [max_sample{k,c}; xaxis(maxind)];
        end
    end
end

%------------ this is full copy of the function above..
function [xaxis,curves,curves_std,sim_curves,sim_curves_std,fraction_sample,cutoff] = get_SAA_samples(handles)

UC = get(handles.C1,'String');

fraction_sample = cell(numel(UC),2); % 2 channels
curves = cell(numel(UC),2);
curves_std = cell(numel(UC),2);
sim_curves = cell(numel(UC),2);
sim_curves_std = cell(numel(UC),2);
cutoff = cell(numel(UC),1);

curves_acc = cell(numel(UC),2);
sim_curves_acc = cell(numel(UC),2);

xaxis = handles.data{1}.SAAstruct{1}.xAxisLabels;

% accumulate
for k=1:numel(handles.data)
        cond_index = find(ismember(UC,handles.conditions{k}));
 
        for m=1:numel(handles.data{k}.SAAstruct)
            if ~isempty(handles.data{k}.SAAstruct{m})
                curve1 = handles.data{k}.SAAstruct{m}.Im1Graph';
                curve2 = handles.data{k}.SAAstruct{m}.Im2Graph';
                sim_curve1 = mean(handles.data{k}.SAAstruct{m}.Im1RandGraph,2);
                sim_curve2 = mean(handles.data{k}.SAAstruct{m}.Im2RandGraph,2);
                
                curves_acc{cond_index,1} = [curves_acc{cond_index,1} curve1];
                curves_acc{cond_index,2} = [curves_acc{cond_index,2} curve2];
                sim_curves_acc{cond_index,1} = [sim_curves_acc{cond_index,1} sim_curve1];
                sim_curves_acc{cond_index,2} = [sim_curves_acc{cond_index,2} sim_curve2];
                %
                f1 = handles.data{k}.SAAstruct{m}.eBarData(1);
                f2 = handles.data{k}.SAAstruct{m}.eBarData(2);
                fraction_sample{cond_index,1} = [fraction_sample{cond_index,1}; f1];
                fraction_sample{cond_index,2} = [fraction_sample{cond_index,2}; f2];
                %        
                cutoff{cond_index} = [cutoff{cond_index}; handles.data{k}.SAAstruct{m}.SAAcutoff];
            end
        end
end

for k=1:numel(UC) % over conditions
        curves{k,1} = median(curves_acc{k,1},2);
        curves{k,2} = median(curves_acc{k,2},2);
        curves_std{k,1} = std(curves_acc{k,1},0,2);
        curves_std{k,2} = std(curves_acc{k,2},0,2);    
            sim_curves{k,1} = median(sim_curves_acc{k,1},2);
            sim_curves{k,2} = median(sim_curves_acc{k,2},2);
            sim_curves_std{k,1} = std(sim_curves_acc{k,1},0,2);
            sim_curves_std{k,2} = std(sim_curves_acc{k,2},0,2);    
        
end


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
