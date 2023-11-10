function varargout = gr_fitting_controls_new_1c(varargin)
% GR_FITTING_CONTROLS_NEW_1C MATLAB code for gr_fitting_controls_new_1c.fig
%      GR_FITTING_CONTROLS_NEW_1C, by itself, creates a new GR_FITTING_CONTROLS_NEW_1C or raises the existing
%      singleton*.
%
%      H = GR_FITTING_CONTROLS_NEW_1C returns the handle to a new GR_FITTING_CONTROLS_NEW_1C or the handle to
%      the existing singleton*.
%
%      GR_FITTING_CONTROLS_NEW_1C('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GR_FITTING_CONTROLS_NEW_1C.M with the given input arguments.
%
%      GR_FITTING_CONTROLS_NEW_1C('Property','Value',...) creates a new GR_FITTING_CONTROLS_NEW_1C or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gr_fitting_controls_new_1c_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gr_fitting_controls_new_1c_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gr_fitting_controls_new_1c

% Last Modified by GUIDE v2.5 09-Nov-2023 16:16:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gr_fitting_controls_new_1c_OpeningFcn, ...
                   'gui_OutputFcn',  @gr_fitting_controls_new_1c_OutputFcn, ...
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


% --- Executes just before gr_fitting_controls_new_1c is made visible.
function gr_fitting_controls_new_1c_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gr_fitting_controls_new_1c (see VARARGIN)

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gr_fitting_controls_new_2c_5p (see VARARGIN)

handles.SMLM_Studio = varargin{1,1};

% Choose default command line output for gr_fitting_controls
%handles.output = hObject;

handles.cutoff = 300; % nm
set(handles.distance_cutoff,'String',num2str(handles.cutoff));
set(handles.gr_shape,'String',{'Gaussian','exponential'});

        handles.colors = zeros(10,3);
        handles.colors(1,:) = [0 0.4470 0.7410];
        handles.colors(2,:) = [0.8500 0.3250 0.0980];
        handles.colors(3,:) = [0.9290 0.6940 0.1250];
        handles.colors(4,:) = [0.4940 0.1840 0.5560];
        handles.colors(5,:) = [0.4660 0.6740 0.1880];
        handles.colors(6,:) = [0.3010 0.7450 0.9330];
        handles.colors(7,:) = [0.6350 0.0780 0.1840];
            handles.colors(8,:) = [0 0 1];
            handles.colors(9,:) = [0 0 0];
            handles.colors(10,:) = [1 0 0];        
        handles.markers = {'o','s','^','d','p','*','+','o','s','^'};
        handles.styles = {'-',':','--','-.','-',':','--','-',':','--'};

handles.ATTR = [];
handles.GR = [];
handles.FITDATA = [];
handles.N_LOCS = [];
handles.AREA = [];

% Choose default command line output for gr_fitting_controls3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if ~isfield(handles.SMLM_Studio,'obj_SMLMdata'), return, end

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
        fov = w; % this is likely wrong, but held for first data reduction cases
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
                    handles.ATTR = [handles.ATTR; [plate condition_index w channel k]];
                    cur_n_locs = cur_data.N_locs;
                    if ~isnumeric(cur_n_locs)                        
                        cur_n_locs = str2double(cur_n_locs);
                    end                    
                    handles.N_LOCS = [handles.N_LOCS; cur_n_locs];                   
                    cur_area = cur_data.Area;
                    if ~isnumeric(cur_area)                        
                        cur_area = str2double(cur_area);
                    end
                    handles.AREA = [handles.AREA; cur_area];                    
                end
            end                                
        end
    end
end
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = gr_fitting_controls_new_1c_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in gr_shape.
function gr_shape_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function gr_shape_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

% --- Executes on button press in fit_and_show_plots.
function fit_and_show_plots_Callback(hObject, eventdata, handles)

    s = get(handles.gr_shape,'String');
    mode = s{get(handles.gr_shape,'Value')};

r = handles.SMLM_Studio.obj_SMLMdata.gr_Object_distance;
%
step = 1;
ctoff = fix(handles.cutoff/step); % nm

n_plots = size(handles.GR,1);

[R,C] = MCED(n_plots);

FITDATA = [];

maxgr = max(handles.GR(:));
mingr = min(handles.GR(:));
maxr = max(r(1:ctoff));

g_bank = cell(n_plots,3);

% fit
hw = waitbar(0,'fitting ROIs, please wait');
for k=1:n_plots
    %
    g_exp = handles.GR(k,:)';
    %
    [g_fit, L, N, n, fval] = fit_gr_1c(r(1:ctoff),g_exp(1:ctoff),handles.N_LOCS(k),handles.AREA(k),mode);
    %
    maxg = (max(g_fit) + max(g_exp))/2;
    relerr = fval/maxg;  
    %
    FITDATA = [FITDATA; [handles.N_LOCS(k), handles.AREA(k), L, N, n, maxg, relerr strcmp(mode,'Gaussian')]];
    %        
    g_bank{k,1} = g_exp;
    g_bank{k,2} = g_fit;
    g_bank{k,3} = relerr;
    if ~isempty(hw), waitbar(k/n_plots); drawnow, end      
end
if ~isempty(hw), delete(hw), drawnow; end

% display
h = figure;
figname = ['1 component fitting, cutoff = ' num2str(ctoff) ' nm, mode - ' mode];
set(h,'Name',figname);
%             
for k=1:n_plots                                                    
   
    g_exp = g_bank{k,1};
    g_fit = g_bank{k,2};
    relerr = g_bank{k,3};
    
    condition_index = handles.ATTR(k,2);                                                    
        
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
    %
    attr = handles.ATTR(k,:); 
    PP = handles.SMLM_Studio.Plate(attr(1));
    CC = handles.SMLM_Studio.Condition{attr(2)};    
    WW = handles.SMLM_Studio.Well{attr(3)};
    cc = handles.SMLM_Studio.Channel{attr(4)};    
    OO = handles.SMLM_Studio.Object{attr(5)};
    sep = ' : ';
    %L1 = [PP sep CC sep WW sep cc sep OO];
    L1 = [CC sep WW sep OO];
    L2 = [num2str(relerr) sep, ... 
         num2str(fix(FITDATA(k,3))) sep, ...
         num2str(fix(FITDATA(k,4)))];      
         legend(ax,{L1,L2},'fontsize',8);       
         axis(ax,[0 maxr, mingr maxgr]);
end
% 
% draw scatterplots
h=figure;
set(h,'Name',['cutoff = ' num2str(ctoff) ' nm, mode ' mode]);
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,2);
ax3 = subplot(2,2,3);
ax4 = subplot(2,2,4);
%
n_conditions = length(handles.SMLM_Studio.Condition);
lwh = 2;
LEGEND = cell(0);
for condition_index=1:n_conditions
    mask = condition_index == handles.ATTR(:,2);
    cd = FITDATA(mask,:); % current data
    %
    %[handles.N_LOCS(k), handles.AREA(k), L, N, n, maxg, relerr]
    Area = cd(:,2)/1e6; % square microns!
    density = cd(:,1)./cd(:,2)*1e6;
    L = cd(:,3);
    N = cd(:,4);
    n = cd(:,5);
    maxg = cd(:,6);
    relerr = cd(:,7);
        
    loglog(ax1,L,N,'color',handles.colors(condition_index,:),'marker','s','linestyle','none','markersize',8,'linewidth',lwh);
    hold(ax1,'on');
    LEGEND = [LEGEND [handles.SMLM_Studio.Condition{condition_index}]];
    %
    semilogy(ax2,N,n,'color',handles.colors(condition_index,:),'marker','s','linestyle','none','markersize',8,'linewidth',lwh);
    hold(ax2,'on');
    %
    semilogx(ax3,Area,density,'color',handles.colors(condition_index,:),'marker','.','linestyle','none','markersize',16,'linewidth',lwh);
    hold(ax3,'on');                       
    
    semilogx(ax4,maxg,relerr,'color',handles.colors(condition_index,:),'marker','.','linestyle','none','markersize',16,'linewidth',lwh);
    hold(ax4,'on'); 
end

hold(ax1,'off');
grid(ax1,'on');
xlabel(ax1,'effective cluster radius [nm]');
ylabel(ax1,'number of localizations per cluster');
legend(ax1,LEGEND,'location','northwest');

hold(ax2,'off');
grid(ax2,'on');
xlabel(ax2,'number of localizations per cluster');
ylabel(ax2,'number of clusters');

hold(ax3,'off');
grid(ax3,'on');
xlabel(ax3,'ROI area [\mu^2]');
ylabel(ax3,'localizations density [1/\mu^2]');

hold(ax4,'off');
grid(ax4,'on');
xlabel(ax4,'max(g(r))');
ylabel(ax4,'fitting error / max(g(r))');

handles.FITDATA = FITDATA;
guidata(hObject, handles);

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
    fh = ancestor(hObject,'figure');     
    delete(fh);

% --- Executes on button press in generate_CSV.
function generate_CSV_Callback(hObject, eventdata, handles)
 %
    if isempty(handles.FITDATA), return, end
    %
    D = [];
    %
       for k = 1:size(handles.FITDATA,1)
            cd = handles.FITDATA(k,:);
            attr = handles.ATTR(k,:); 
            PP = {handles.SMLM_Studio.Plate(attr(1))};
            CC = {handles.SMLM_Studio.Condition{attr(2)}};    
            WW = {handles.SMLM_Studio.Well{attr(3)}};
            cc = {handles.SMLM_Studio.Channel{attr(4)}};    
            OO = {handles.SMLM_Studio.Object{attr(5)}};

            rec = [PP CC WW cc OO num2cell(cd)];
            D = [D; rec];
       end 
       D
       %
       i1 = cell2mat(D(:,13));
       type1=cell(size(i1));
       for k=1:length(i1)
            if i1(k), type1{k}='Gauss'; else type1{k}='exp'; end
       end
       %
       D(:,13) = type1;
       %         
       caption = {'plate','condition','well','channel','object', ...
       '#locs','Area','Z','N','n','max(g)','fit_err/max(g)','F'};    

D = [caption; D];

xlstempname = [tempname '.csv'];
cell2csv(xlstempname,D);
if ispc 
    winopen(xlstempname);
else
    open(xlstempname);
end    

%-----------------------------------------------------------------
function distance_cutoff_Callback(hObject, eventdata, handles)
max_r = max(handles.SMLM_Studio.obj_SMLMdata.gr_Object_distance);
%
value = fix(str2double(get(hObject,'String')));
if ~isnan(value) && value>=40 && value <= max_r
    handles.cutoff = value;
    guidata(hObject,handles);
else
    value = handles.cutoff;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% Hints: get(hObject,'String') returns contents of distance_cutoff as text
%        str2double(get(hObject,'String')) returns contents of distance_cutoff as a double


% --- Executes during object creation, after setting all properties.
function distance_cutoff_CreateFcn(hObject, eventdata, handles)
max_r = max(handles.SMLM_Studio.obj_SMLMdata.gr_Object_distance);
%
value = fix(str2double(get(hObject,'String')));
if ~isnan(value) && value>=40 && value <= max_r
    handles.cutoff = value;
    guidata(hObject,handles);
else
    value = handles.cutoff;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

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