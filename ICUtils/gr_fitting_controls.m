function varargout = gr_fitting_controls(varargin)
% GR_FITTING_CONTROLS MATLAB code for gr_fitting_controls.fig
%      GR_FITTING_CONTROLS, by itself, creates a new GR_FITTING_CONTROLS or raises the existing
%      singleton*.
%
%      H = GR_FITTING_CONTROLS returns the handle to a new GR_FITTING_CONTROLS or the handle to
%      the existing singleton*.
%
%      GR_FITTING_CONTROLS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GR_FITTING_CONTROLS.M with the given input arguments.
%
%      GR_FITTING_CONTROLS('Property','Value',...) creates a new GR_FITTING_CONTROLS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gr_fitting_controls_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gr_fitting_controls_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gr_fitting_controls

% Last Modified by GUIDE v2.5 01-Nov-2022 13:47:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gr_fitting_controls_OpeningFcn, ...
                   'gui_OutputFcn',  @gr_fitting_controls_OutputFcn, ...
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


% --- Executes just before gr_fitting_controls is made visible.
function gr_fitting_controls_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gr_fitting_controls (see VARARGIN)

handles.SMLM_Studio = varargin{1,1};

if ~isfield(handles.SMLM_Studio,'obj_SMLMdata'), return, end

% Choose default command line output for gr_fitting_controls
%handles.output = hObject;

handles.cutoff = 300; % nm
set(handles.distance_cutoff,'String',num2str(handles.cutoff));
set(handles.number_of_parameters,'String',{'4','5'});
set(handles.shape_primitives,'String',{'Gaussian','exponential','mixed'});

        handles.colors = zeros(7,3);
        handles.colors(1,:) = [0 0.4470 0.7410];
        handles.colors(2,:) = [0.8500 0.3250 0.0980];
        handles.colors(3,:) = [0.9290 0.6940 0.1250];
        handles.colors(4,:) = [0.4940 0.1840 0.5560];
        handles.colors(5,:) = [0.4660 0.6740 0.1880];
        handles.colors(6,:) = [0.3010 0.7450 0.9330];
        handles.colors(7,:) = [0.6350 0.0780 0.1840];
        handles.markers = {'o','s','^','d','p','*','+'};
        handles.styles = {'-',':','--','-.','-',':','--'};

handles.ATTR = [];
handles.GR = [];
handles.FIT = [];
handles.N_LOCS = [];
handles.AREA = [];
        
d = handles.SMLM_Studio.obj_SMLMdata;
data = d.data;
[n_plates, ...
n_conditions, ...
n_wells, ...
n_fovs, ...
n_channels, ...
n_max_objects] = size(data);  
%
for plate = 1:n_plates
    for w=1:n_wells        
        well_is_empty = true;
        try
            well_is_empty = isempty(handles.SMLM_Studio.well_to_condition_map(num2str(w)));
        catch
            well_is_empty = isempty(handles.SMLM_Studio.well_to_condition_map(handles.SMLM_Studio.Well{w}));
        end
        %
        fov = w; % this is lkely wrong, but held for first data reduction cases
        %
        for k = 1:n_max_objects
            for channel = 1:n_channels
                %
                try
                    condition_index = handles.SMLM_Studio.well_to_condition_index_map(num2str(w));
                catch
                    condition_index = handles.SMLM_Studio.well_to_condition_index_map(handles.SMLM_Studio.Well{w});
                end
                %
                if ~isempty(condition_index) && 0 ~= condition_index
                    cur_data = data{plate,condition_index,w,fov,channel,k};
                    if isempty(cur_data), continue, end
                    disp([plate w fov k channel])                    
                    %
                    g_exp = cur_data.gr;
                    handles.GR = [handles.GR; g_exp'];
                    handles.ATTR = [handles.ATTR; [plate w fov k channel]];
                    handles.N_LOCS = [handles.N_LOCS; cur_data.N_locs];
                    handles.AREA = [handles.AREA; cur_data.Area];                    
                end
            end                                
        end
    end
end        
                                
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gr_fitting_controls wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gr_fitting_controls_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


% --- Executes on selection change in distance_cutoff.
function distance_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to distance_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns distance_cutoff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from distance_cutoff

max_r = max(handles.SMLM_Studio.obj_SMLMdata.gr_Object_distance);
%
value = fix(str2double(get(hObject,'String')));
if ~isnan(value) && value>50 && value <= max_r
    handles.cutoff = value;
    guidata(hObject,handles);
else
    value = handles.cutoff;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function distance_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in number_of_parameters.
function number_of_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns number_of_parameters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from number_of_parameters


% --- Executes during object creation, after setting all properties.
function number_of_parameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in shape_primitives.
function shape_primitives_Callback(hObject, eventdata, handles)
% hObject    handle to shape_primitives (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns shape_primitives contents as cell array
%        contents{get(hObject,'Value')} returns selected item from shape_primitives


% --- Executes during object creation, after setting all properties.
function shape_primitives_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shape_primitives (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fit_and_show_plots.
function fit_and_show_plots_Callback(hObject, eventdata, handles)
% hObject    handle to fit_and_show_plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
    s = get(handles.number_of_parameters,'String');
num_param = fix(str2double(s{get(handles.number_of_parameters,'Value')}));
    s = get(handles.shape_primitives,'String');
mode = s{get(handles.shape_primitives,'Value')};

r = handles.SMLM_Studio.obj_SMLMdata.gr_Object_distance;
%
step = 1;
ctoff = fix(handles.cutoff/step); % nm

n_plots = size(handles.GR,1);

[R,C] = MCED(n_plots);

figure('units','normalized','outerposition',[0 0 1 1],'name','g(r) fittings');
for k=1:n_plots
    k
    g_exp = handles.GR(k,:)';
    condition_index = handles.ATTR(k,2);
    [g_fit, L1, L2, n1, n2, p1, p2, N1, N2, fval] = fit_gr(r(1:ctoff),g_exp(1:ctoff),...
        handles.N_LOCS(k),handles.AREA(k),mode,num_param);
    %    
    ax = subplot(R,C,k);
                                        
    semilogy(ax,r(1:ctoff),g_exp(1:ctoff),'color',handles.colors(condition_index,:),'marker',handles.markers{condition_index});
    hold(ax,'on');
    semilogy(ax,r(1:ctoff),g_fit,'k:','linewidth',2);
    hold(ax,'on');                        

    hold(ax,'off');
    grid(ax,'on');
    if k==n_plots
        xlabel(ax,'distance [nm]');
    end
    if k==1
        ylabel(ax,'g(r)');
    end
    legend(ax,[handles.SMLM_Studio.Condition{condition_index} ' : ' num2str(fval)]);
end

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
    fh = ancestor(hObject,'figure');     
    delete(fh);

% --- Executes on button press in fit_and_generate_xls.
function fit_and_generate_xls_Callback(hObject, eventdata, handles)
% hObject    handle to fit_and_generate_xls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% -------------------------------------------------------------
function [med1,med2] = MCED(N_)
%MCED maximum close to each other divisors
%   
N = abs(round(N_));
med1 = [];
med2 = [];

H = floor(sqrt(N));

for k=H:-1:1
    if 0==rem(N,k)
       med1 = k;
       med2 = N/k;
       break;
    end
end
