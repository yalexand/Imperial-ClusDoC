
classdef ClusDoC_analysis_controller < handle 
    
    % 2020 Imperial College London.
    % this module is an adaptation of ClusDoC software
   
    properties(Constant)
        data_settings_filename = 'ClusDoC_settings.xml';
    end
    
    properties(SetObservable = true)
	
    pixelSizenm = 150;
    WindSTORM_Sigmapix = 1;

    % Default RipleyK settings
    RipleyK = struct( ...
        'Start',0, ...
        'End',1000, ...
        'Step',10, ...
        'MaxSampledPts', 1e4);
    
    % Default DBSCAN parameters
    DBSCAN = struct( ...
        'epsilon',20, ...
        'minPts',3, ...
        'UseLr_rThresh',true, ...
        'Lr_rThreshRad',20, ...
        'SmoothingRad',15, ...
        'Cutoff',10, ...
        'threads',12, ...
        'DoStats',true, ...
        'Outputfolder',['C:' filesep]);
 
        Square_ROIs_Auto_anm = 1400;
        Square_ROIs_Auto_qthresh = [.5 .5];
        Square_ROIs_Auto_maxNrois = 20; 
        Square_ROIs_Auto_method = 'channel';
        Square_ROIs_Auto_LNT = [200 200]; % loc number thershold - may be interpreted as per_micron or absolute
        Square_ROIs_Auto_HNT = [20000 20000];
        Square_ROIs_Auto_verbose = false;
        
        Align_channels_nmppix = 4;
        Align_channels_method = 'Matlab_multimodal';
    
    % Initialize structure to pass values between GUI components
    CellData = {2,1}; % better it rename this as FOVData
    
    % original data
    Original_from_file_data = {2,1};
    Original_from_file_header = {2,1};
    
    pathName = cell(2,1);
    fileName = cell(2,1);
    
    NDataColumns = [];
    
    SizeX = 256; % # original image pixels as on acquisition
    SizeY = 256;

    Nchannels = 1;
    
    ROICoordinates = {};
    
    Outputfolder = [];
    
    ClusterTable = [];
    
    Chan1Color = [1 0 0]; %red   
    
    sgm = []; % segmentation image
    sgm_base = [];
    
    roi_to_object_lut = []; % square ROIs are chosen within "object" 
    
    save_DBSCAN_clusters_images = false; % to reduce output results storage size
    
%%%%%%%%%%%%%%%%%%%%%%%% DoC                     
    % DoC DBSCAN parameters
    DoC_dbscanParams_ch1 = struct( ...
        'DoCThreshold',0.4, ...
        'epsilon',20, ...
        'minPts',3, ...
        'UseLr_rThresh',true, ...
        'Lr_rThreshRad',20, ...
        'SmoothingRad',15, ...
        'Cutoff',10, ...
        'threads',12, ...
        'DoStats',true, ...
        'Outputfolder',['C:' filesep]);                  
    DoC_dbscanParams_ch2 = struct( ...
        'DoCThreshold',0.4, ...
        'epsilon',20, ...
        'minPts',3, ...
        'UseLr_rThresh',true, ...
        'Lr_rThreshRad',20, ...
        'SmoothingRad',15, ...
        'Cutoff',10, ...
        'threads',12, ...
        'DoStats',true, ...
        'Outputfolder',['C:' filesep]);           
    
        DoC_Lr_rRad = 20;
        DoC_Rmax = 500;
        DoC_Step = 10;
        DoC_ColoThres = 0.4;            
        DoC_NbThresh = 10;
%%%%%%%%%%%%%%%%%%%%%%%% DoC

        %
        GR_ASSISTED_CLUSTERING = struct( ...
        'PERFORM',false, ...
        'r_start',750, ...
        'dr',40, ...
        'MODE','2-comp Gauss', ...
        'minPts1',5, ...
        'minPts2',5, ...        
        'min_Z2_Z1_ratio',2, ...
        'max_Z2_allowed',500, ...
        'DBSCAN_epsilon',40, ...
        'DBSCAN_SmoothingRad',40, ...
        'Small_Clusters_Size_Tolerance_Factor',0.75, ...
        'imopen_Z1_excess_factor',2.5, ...
        'sampling_Z1_excess_factor',1.2, ...
        'sampling_Z2_excess_factor',1.2);
        %

        Plate = [];
        Well = [];
        Condition = [];
        
    end                    
    
    properties(Transient)      
        DefaultDirectory = ['C:' filesep];
        RootDirectory = ['C:' filesep];
      end    
        
    properties(Transient,Hidden)    
    end
    
    events    
    end
            
    methods
        
        function obj = ClusDoC_analysis_controller(varargin)
            this_dir = pwd;
            if ~isdeployed
                resolve_xlwrite_path_issue();
            end
            obj.load_settings([this_dir filesep obj.data_settings_filename]);
        end

%-------------------------------------------------------------------------%
function load_settings(obj,fname,~)
    try
             if exist(fname,'file') 
                [ settings, ~ ] = xml_read (fname);
                    obj.pixelSizenm = settings.pixelSizenm;
                    obj.WindSTORM_Sigmapix = settings.WindSTORM_Sigmapix;

                    obj.RipleyK.Start = settings.RipleyK_Start;
                    obj.RipleyK.End = settings.RipleyK_End;
                    obj.RipleyK.Step = settings.RipleyK_Step;
                    obj.RipleyK.MaxSampledPts = settings.RipleyK_MaxSampledPts;

                    obj.DBSCAN.epsilon = settings.DBSCAN_epsilon;
                    obj.DBSCAN.minPts = settings.DBSCAN_minPts;
                    obj.DBSCAN.UseLr_rThresh =settings.DBSCAN_UseLr_rThresh;
                    obj.DBSCAN.Lr_rThreshRad = settings.DBSCAN_Lr_rThreshRad;
                    obj.DBSCAN.SmoothingRad = settings.DBSCAN_SmoothingRad;
                    obj.DBSCAN.Cutoff = settings.DBSCAN_Cutoff;
                    obj.DBSCAN.threads = settings.DBSCAN_threads;
                    obj.DBSCAN.DoStats = settings.DBSCAN_DoStats;
                    
                    obj.Square_ROIs_Auto_anm = settings.Square_ROIs_Auto_anm; % square side
                    obj.Square_ROIs_Auto_qthresh = settings.Square_ROIs_Auto_qthresh; % brightness threshold
                    obj.Square_ROIs_Auto_maxNrois = settings.Square_ROIs_Auto_maxNrois;
                    obj.Square_ROIs_Auto_method = settings.Square_ROIs_Auto_method;
                    obj.Square_ROIs_Auto_LNT = settings.Square_ROIs_Auto_LNT;
                    obj.Square_ROIs_Auto_HNT = settings.Square_ROIs_Auto_HNT;
                    
                    obj.SizeX = settings.SizeX;
                    obj.SizeY = settings.SizeY;                           
        
                    obj.Align_channels_nmppix = settings.Align_channels_nmppix;
                    obj.Align_channels_method = settings.Align_channels_method;
                    
%%%%%%%%%%%%%%%%%%%%%%%% DoC                     
    % DoC DBSCAN parameters
                    obj.DoC_dbscanParams_ch1.DoCThreshold = settings.DoC_dbscanParams_ch1_DoCThreshold;
                    obj.DoC_dbscanParams_ch1.epsilon = settings.DoC_dbscanParams_ch1_epsilon;
                    obj.DoC_dbscanParams_ch1.minPts = settings.DoC_dbscanParams_ch1_minPts;
                    obj.DoC_dbscanParams_ch1.UseLr_rThresh = settings.DoC_dbscanParams_ch1_UseLr_rThresh;
                    obj.DoC_dbscanParams_ch1.Lr_rThreshRad = settings.DoC_dbscanParams_ch1_Lr_rThreshRad;
                    obj.DoC_dbscanParams_ch1.SmoothingRad = settings.DoC_dbscanParams_ch1_SmoothingRad;
                    obj.DoC_dbscanParams_ch1.Cutoff = settings.DoC_dbscanParams_ch1_Cutoff;
                    obj.DoC_dbscanParams_ch1.threads = settings.DoC_dbscanParams_ch1_threads;
                    obj.DoC_dbscanParams_ch1.DoStats = settings.DoC_dbscanParams_ch1_DoStats;
                    obj.DoC_dbscanParams_ch1.Outputfolder = settings.DoC_dbscanParams_ch1_Outputfolder;

                    obj.DoC_dbscanParams_ch2.DoCThreshold = settings.DoC_dbscanParams_ch2_DoCThreshold;
                    obj.DoC_dbscanParams_ch2.epsilon = settings.DoC_dbscanParams_ch2_epsilon;
                    obj.DoC_dbscanParams_ch2.minPts = settings.DoC_dbscanParams_ch2_minPts;
                    obj.DoC_dbscanParams_ch2.UseLr_rThresh = settings.DoC_dbscanParams_ch2_UseLr_rThresh;
                    obj.DoC_dbscanParams_ch2.Lr_rThreshRad = settings.DoC_dbscanParams_ch2_Lr_rThreshRad;
                    obj.DoC_dbscanParams_ch2.SmoothingRad = settings.DoC_dbscanParams_ch2_SmoothingRad;
                    obj.DoC_dbscanParams_ch2.Cutoff = settings.DoC_dbscanParams_ch2_Cutoff;
                    obj.DoC_dbscanParams_ch2.threads = settings.DoC_dbscanParams_ch2_threads;
                    obj.DoC_dbscanParams_ch2.DoStats = settings.DoC_dbscanParams_ch2_DoStats;
                    obj.DoC_dbscanParams_ch2.Outputfolder = settings.DoC_dbscanParams_ch2_Outputfolder;

                    obj.DoC_Lr_rRad = settings.DoC_Lr_rRad;
                    obj.DoC_Rmax = settings.DoC_Rmax;
                    obj.DoC_Step = settings.DoC_Step;
                    obj.DoC_ColoThres = settings.DoC_ColoThres;   
                    obj.DoC_NbThresh = settings.DoC_NbThresh;
%%%%%%%%%%%%%%%%%%%%%%%% DoC                    
                                       
             end
    catch err
        disp('cannot load settings file, exiting..');
        disp(err.message);
    end
end
%-------------------------------------------------------------------------%
function save_settings(obj,fname,~)
    try
                    settings = [];
                    settings.pixelSizenm = obj.pixelSizenm;
                    settings.WindSTORM_Sigmapix = obj.WindSTORM_Sigmapix;

                    settings.RipleyK_Start = obj.RipleyK.Start;
                    settings.RipleyK_End = obj.RipleyK.End;
                    settings.RipleyK_Step = obj.RipleyK.Step;
                    settings.RipleyK_MaxSampledPts = obj.RipleyK.MaxSampledPts;

                    settings.DBSCAN_epsilon = obj.DBSCAN.epsilon;
                    settings.DBSCAN_minPts = obj.DBSCAN.minPts;
                    settings.DBSCAN_UseLr_rThresh = obj.DBSCAN.UseLr_rThresh;
                    settings.DBSCAN_Lr_rThreshRad = obj.DBSCAN.Lr_rThreshRad;
                    settings.DBSCAN_SmoothingRad = obj.DBSCAN.SmoothingRad;
                    settings.DBSCAN_Cutoff = obj.DBSCAN.Cutoff;
                    settings.DBSCAN_threads = obj.DBSCAN.threads;
                    settings.DBSCAN_DoStats = obj.DBSCAN.DoStats;
                    
                    settings.Square_ROIs_Auto_anm = obj.Square_ROIs_Auto_anm; % square side
                    settings.Square_ROIs_Auto_qthresh = obj.Square_ROIs_Auto_qthresh; % brightness threshold
                    settings.Square_ROIs_Auto_maxNrois = obj.Square_ROIs_Auto_maxNrois;
                    settings.Square_ROIs_Auto_method = obj.Square_ROIs_Auto_method;
                    settings.Square_ROIs_Auto_LNT = obj.Square_ROIs_Auto_LNT;
                    settings.Square_ROIs_Auto_HNT = obj.Square_ROIs_Auto_HNT;
                    
                    settings.SizeX = obj.SizeX;
                    settings.SizeY = obj.SizeY; 
                    
                    settings.Align_channels_nmppix = obj.Align_channels_nmppix;
                    settings.Align_channels_method = obj.Align_channels_method; 
                    
%%%%%%%%%%%%%%%%%%%%%%%% DoC                     
    % DoC DBSCAN parameters
                    settings.DoC_dbscanParams_ch1_DoCThreshold = obj.DoC_dbscanParams_ch1.DoCThreshold;
                    settings.DoC_dbscanParams_ch1_epsilon = obj.DoC_dbscanParams_ch1.epsilon;
                    settings.DoC_dbscanParams_ch1_minPts = obj.DoC_dbscanParams_ch1.minPts;
                    settings.DoC_dbscanParams_ch1_UseLr_rThresh = obj.DoC_dbscanParams_ch1.UseLr_rThresh;
                    settings.DoC_dbscanParams_ch1_Lr_rThreshRad = obj.DoC_dbscanParams_ch1.Lr_rThreshRad;
                    settings.DoC_dbscanParams_ch1_SmoothingRad = obj.DoC_dbscanParams_ch1.SmoothingRad;
                    settings.DoC_dbscanParams_ch1_Cutoff = obj.DoC_dbscanParams_ch1.Cutoff;
                    settings.DoC_dbscanParams_ch1_threads = obj.DoC_dbscanParams_ch1.threads;
                    settings.DoC_dbscanParams_ch1_DoStats = obj.DoC_dbscanParams_ch1.DoStats;
                    settings.DoC_dbscanParams_ch1_Outputfolder = obj.DoC_dbscanParams_ch1.Outputfolder;

                    settings.DoC_dbscanParams_ch2_DoCThreshold = obj.DoC_dbscanParams_ch2.DoCThreshold;
                    settings.DoC_dbscanParams_ch2_epsilon = obj.DoC_dbscanParams_ch2.epsilon;
                    settings.DoC_dbscanParams_ch2_minPts = obj.DoC_dbscanParams_ch2.minPts;
                    settings.DoC_dbscanParams_ch2_UseLr_rThresh = obj.DoC_dbscanParams_ch2.UseLr_rThresh;
                    settings.DoC_dbscanParams_ch2_Lr_rThreshRad = obj.DoC_dbscanParams_ch2.Lr_rThreshRad;
                    settings.DoC_dbscanParams_ch2_SmoothingRad = obj.DoC_dbscanParams_ch2.SmoothingRad;
                    settings.DoC_dbscanParams_ch2_Cutoff = obj.DoC_dbscanParams_ch2.Cutoff;
                    settings.DoC_dbscanParams_ch2_threads = obj.DoC_dbscanParams_ch2.threads;
                    settings.DoC_dbscanParams_ch2_DoStats = obj.DoC_dbscanParams_ch2.DoStats;
                    settings.DoC_dbscanParams_ch2_Outputfolder = obj.DoC_dbscanParams_ch2.Outputfolder;

                    settings.DoC_Lr_rRad = obj.DoC_Lr_rRad;
                    settings.DoC_Rmax = obj.DoC_Rmax;
                    settings.DoC_Step = obj.DoC_Step;
                    settings.DoC_ColoThres = obj.DoC_ColoThres;   
                    settings.DoC_NbThresh = obj.DoC_NbThresh;
%%%%%%%%%%%%%%%%%%%%%%%% DoC               
                                                                               
                    xml_write(fname,settings);
    catch err
        disp('cannot save settings file');
        disp(err.message);
    end
end
%-------------------------------------------------------------------------%
        function delete(obj)           
        end       
%-------------------------------------------------------------------------%                
        function clear_all(obj,~)
        end
%-------------------------------------------------------------------------%
        function Load_Data(obj,fileName,pathName,chan,~)
            
            obj.fileName{chan} = fileName;
            obj.pathName{chan} = pathName;
        
            goodZENFile = checkZenFile(fullfile(obj.pathName{chan},obj.fileName{chan}));
      
            if goodZENFile || contains(obj.fileName{chan},'.csv') || contains(obj.fileName{chan},'.txt')
                if goodZENFile || contains(obj.fileName{chan},'.csv')
                    [importData,original_header,original_data] = ImportFile(fullfile(obj.pathName{chan},obj.fileName{chan}),obj);
                else % QC_STORM
                    [importData,original_header,original_data] = Import_QC_STORM_File(fullfile(obj.pathName{chan},obj.fileName{chan}),obj);   
                end
                                
                obj.CellData{chan} = [importData zeros(size(importData, 1), 8)];

%                 obj.CellData{k}(any(isnan(obj.CellData{k}), 2), :) = []; % protection against incomplete line writing in ZEN export
                                                                                 % This breaks import for ThunderSTORMConcatenator output                                                                                  % Commenting this out to allow that format to work.
                obj.NDataColumns = size(importData, 2);
                obj.CellData{chan}(:,obj.NDataColumns + 2) = 1; % All data is in mask until set otherwise
                                             
%                 obj.CellData{chan}(any(obj.CellData{chan}(:, 5) > obj.SizeX), : )= [];
%                 obj.CellData{chan}(any(obj.CellData{chan}(:, 6) > obj.SizeY), : )= [];                
%                 obj.CellData{chan}(any(obj.CellData{chan}(:, 5:6) < 0), : )= [];
                
                obj.Nchannels = max(chan,obj.Nchannels); %min([numel(unique(obj.CellData{chan}(:,12))), 2]); % cap import to 2 channels ever
                
                obj.Original_from_file_data{chan} = original_data;
                obj.Original_from_file_header{chan} = original_header;
                
            else
                
                fprintf(1, 'File not in accepted coordinate table format.\nSkipping %s\n', fullfile(obj.pathName{chan},obj.fileName{chan}));

            end            
            
                %-------------------------------------------------------------------------%
                % Checking function for selected files
                function isGood = checkZenFile(fName)

                    fID = fopen(fName, 'r');
                    firstLine = fgetl(fID);
                    nTabs = length(strfind(firstLine, sprintf('\t')));
                    firstEntry = firstLine(1:5);
                    fclose(fID);

                    if ismember(nTabs, [11, 12, 13, 14]) && strcmp(firstEntry, 'Index')
                        isGood = true;
                    elseif ismember(nTabs, [23, 25]) && strcmp(firstEntry, 'Chann')
                        isGood = true;
                        % Is good Nikon file, which will get interpreted into Zeiss
                        % format in Import1File
                    elseif ismember(nTabs, 0) 

                            fID = fopen(fName, 'r');
                            firstLine = fgetl(fID);
                            nTabs = length(strfind(firstLine, sprintf(',')));
                            firstEntry = firstLine(1:5);
                            fclose(fID);

                            if ismember(nTabs, [9, 10]) && strcmp(firstEntry, '"id",')
                                isGood = true;
                                % Is good ThunderSTORM file
                            else
                                isGood = false;
                            end
                    else
                        isGood = false;
                    end
                end % checkZenFile
                %             
        end
%-------------------------------------------------------------------------%
    function Analyze_ROIs_DBSCAN(obj,verbose,~)
  
        try
            
        for chan = 1 : obj.Nchannels
               
               fname = strrep(obj.fileName{chan},'.csv','');
               fname = strrep(fname,'.txt','');
               DBSCAN_out_dirname = [obj.Outputfolder filesep fname '_channel_' num2str(chan) '_ClusDoC_Results' filesep 'DBSCAN'];         
               if ~exist(DBSCAN_out_dirname,'dir')
                   mkdir( fullfile(obj.Outputfolder,[fname '_channel_' num2str(chan) '_ClusDoC_Results'],'DBSCAN'));
                   mkdir(DBSCAN_out_dirname,'Cluster_maps');
                   mkdir(DBSCAN_out_dirname,'Cluster_density_maps');
               end
                          
        obj.DBSCAN.Outputfolder = DBSCAN_out_dirname;  
               
        cellROIPair = [];                  
               
        Result = cell(length(obj.ROICoordinates),1);
        ClusterSmoothTable = cell(length(obj.ROICoordinates),1);
                
                    if verbose
                        %figure;
                        figure('units','normalized','outerposition',[0 0 1 1],'name','ROIs as they go..');
                        ax = gca;
                        plot(ax,obj.CellData{chan}(:,5),obj.CellData{chan}(:,6),'b.');
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                    end
                    
                    for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                    x_roi=roi(:,1);
                    y_roi=roi(:,2);
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                    dataCropped = obj.CellData{chan}(whichPointsInROI,:);

                    if verbose
                        color = [rand rand 0];
                        plot(ax,dataCropped(:,5),dataCropped(:,6),'marker','.','color',color,'linestyle','none');
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                    end
                     
                    c = 1;
                        if ~isempty(dataCropped)
                        
                            % DBSCANHandler(Data, DBSCANParams, varargin)
                            %         p = varargin{1}; % Labeling only
                            %         q = varargin{2}; % Labeling only
                            %         display1 = varargin{3};
                            %         display2 = varargin{4};
                            %         clusterColor = varargin{5}
                            
                            clusterColor = obj.Chan1Color;
%                             [~, ClusterSmoothTable{roiInc}, ~, classOut, ~, ~, ~, Result{roiInc,c}] = ...
%                                 DBSCANHandler_YA(dataCropped(dataCropped(:,12) == chan, 5:6), obj.DBSCAN, c, roiInc, ...
%                                 true, true, clusterColor, dataCropped(dataCropped(:,12) == chan, obj.NDataColumns + 2));
                            try
                                [~, ClusterSmoothTable{roiInc,c}, ~, classOut, ~, ~, ~, Result{roiInc,c}] = ...
                                    obj.DBSCANHandler_YA(dataCropped(dataCropped(:,12) == 1, 5:6), obj.DBSCAN, c, roiInc, ...
                                    true, true, clusterColor, dataCropped(dataCropped(:,12) == 1, obj.NDataColumns + 2));                               
                                
                                %obj.CellData(whichPointsInROI & (obj.CellData(:,12) == chan), obj.NDataColumns + 3) = classOut;
                                obj.CellData{chan}(whichPointsInROI & (obj.CellData{chan}(:,12) == 1), obj.NDataColumns + 3) = classOut;
                                
                                % obj.ClusterTable = AppendToClusterTable(obj.ClusterTable, chan, c, roiInc, ClusterSmoothTable{roiInc, c}, classOut);
                                obj.ClusterTable = AppendToClusterTable(obj.ClusterTable, 1, c, roiInc, ClusterSmoothTable{roiInc, c}, classOut);                                 
                                
                            catch err                        
                                fprintf(1, 'ROI %d - error.  Skipping.\n', roiInc);
                                disp(err.message);
                                ClusterSmoothTable{roiInc,c} = [];
                                classOut = [];
                                Result{roiInc,c} = [];
                            end                           
                                
                        else
                            % Have chosen an empty region as ROI                            
                            fprintf(1, 'Cell %d - ROI %d is empty.  Skipping.\n', 1, roiInc);                           
                            ClusterSmoothTable{roiInc,c} = [];
                            classOut = [];
                            Result{roiInc,c} = [];
                        end
                        
                        cellROIPair = [cellROIPair; c, roiInc, roi(1,1), roi(1,2), polyarea(roi(:,1), roi(:,2))];
                        
                    end % ROI
                    
                    if verbose, hold(ax,'off'), end

                if ~all(cellfun(@isempty, Result))
                    
                   try
                    %ExportDBSCANDataToExcelFiles(cellROIPair, Result, obj.DBSCAN.Outputfolder, chan);
                    ExportDBSCANDataToExcelFiles_YA(cellROIPair, Result, obj.DBSCAN.Outputfolder, 1);
                   catch err
                       disp('error in function ExportDBSCANDataToExcelFiles_YA');
                       disp(err.message);
                   end
                else
                    fprintf(1, 'All cells and ROIs empty.  Skipping export.\n');
                end
                
               save(fullfile(obj.DBSCAN.Outputfolder,'DBSCAN_Cluster_Result.mat'),'ClusterSmoothTable','Result','-v7.3');                
              
        end
               
        catch mErr
                
           disp('DBSCAN processing exited with errors.');
           rethrow(mErr);

        end

                        %-------------------------------------------------------------------------%
                        function clusterTableOut = AppendToClusterTable(clusterTable, Ch, cellIter, roiIter, ClusterCh, classOut)

                            try 
                                if isempty(clusterTable)
                                    oldROIRows = [];
                                else
                                    oldROIRows = (cellIter == clusterTable(:,1)) & (roiIter == clusterTable(:,2)) & (Ch == clusterTable(:,3));
                                end

                                if any(oldROIRows)

                                    % Clear out the rows that were for this ROI done previously
                                    clusterTable(oldROIRows, :) = [];

                                end

                                % Add new data to the clusterTable
                                appendTable = nan(length(ClusterCh), 16);
                                appendTable(:, 1) = cellIter; % CurrentROI
                                appendTable(:, 2) = roiIter; % CurrentROI
                                appendTable(:, 3) = Ch; % Channel

                                appendTable(:, 4) = cellfun(@(x) x.ClusterID, ClusterCh); % ClusterID
                                appendTable(:, 5) = cell2mat(cellfun(@(x) size(x.Points, 1), ClusterCh, 'uniformoutput', false)); % NPoints
                                appendTable(:, 6) = cellfun(@(x) x.Nb, ClusterCh); % Nb

                                if ~isempty(ClusterCh)

                                    if isfield(ClusterCh{1}, 'MeanDoC')
                                        appendTable(:, 7) = cellfun(@(x) x.MeanDoC, ClusterCh); % MeanDoCScore
                                    end

                                    if isfield(ClusterCh{1}, 'AvRelativeDensity')
                                        appendTable(:, 11) = cellfun(@(x) x.AvRelativeDensity, ClusterCh); % AvRelativeDensity
                                        appendTable(:, 12) = cellfun(@(x) x.Mean_Density, ClusterCh); % MeanDensity
                                    end

                                    if isfield(ClusterCh{1}, 'Nb_In')
                                        appendTable(:, 13) = cellfun(@(x) x.Nb_In, ClusterCh); % Nb_In
                                    end

                                end

                                appendTable(:, 8) = cellfun(@(x) x.Area, ClusterCh); % Area
                                appendTable(:, 9) = cellfun(@(x) x.Circularity, ClusterCh); % Circularity
                                appendTable(:, 10) = cellfun(@(x) x.TotalAreaDensity, ClusterCh); % TotalAreaDensity

                                appendTable(:, 14) = cellfun(@(x) x.NInsideMask, ClusterCh); % NPointsInsideMask
                                appendTable(:, 15) = cellfun(@(x) x.NOutsideMask, ClusterCh); % NPointsInsideMask
                                appendTable(:, 16) = cellfun(@(x) x.Elongation, ClusterCh); % Elongation                                

                                clusterTableOut = [clusterTable; appendTable];

                            catch mError
                                assignin('base', 'ClusterCh', ClusterCh);
                        %         assignin('base', 'clusterIDList', clusterIDList);
                                assignin('base', 'appendTable', appendTable);
                                assignin('base', 'classOut', classOut);

                                rethrow(mError);

                            end
                        end
                        %-------------------------------------------------------------------------%        
    end
%-------------------------------------------------------------------------%
    function Analyze_ROIs_RipleyK(obj,~)
  
        try
    
            % Ripley K calculation
            % Iterate through cells + ROIs

            % Code inside loop is from RipleyKmultiData_GUIFunV2.m 
            % Moving here to create more reasonable workflow
            % create the output folder 'RipleyKGUI_Result

               for chan = 1 : obj.Nchannels
                   
                   fname = strrep(obj.fileName{chan},'.csv','');
                   fname = strrep(fname,'.txt','');
                   RipleyK_out_dirname = [obj.Outputfolder filesep fname '_channel_' num2str(chan) '_ClusDoC_Results' filesep  'RipleyK'];
                   if ~exist(RipleyK_out_dirname,'dir')
                       mkdir( fullfile(obj.Outputfolder,[fname '_channel_' num2str(chan) '_ClusDoC_Results'],'RipleyK'));
                       mkdir(RipleyK_out_dirname,'RipleyK_Plots');
                       mkdir(RipleyK_out_dirname,'RipleyK_Results');
                   end
                                      
                   [~] = obj.RipleyKHandler(RipleyK_out_dirname,chan);
            
               end                       
            
        catch mError
            
            disp('Ripley K processing exited with errors');
            %rethrow(mError);
            disp(mError.message);
             
        end        
    end
    %-------------------------------------------------------------------------%   
    function valOut = RipleyKHandler(obj,Fun_OutputFolder_name,chan,~)
        
        if isempty(obj.ROICoordinates), return, end
        
        Start = obj.RipleyK.Start;
        End = obj.RipleyK.End;
        Step = obj.RipleyK.Step;
        MaxSampledPts = obj.RipleyK.MaxSampledPts;  

            % Handler function for RipleyK calculations
                ArrayHeader = [{'r'},{'L(r)-r'}];
                
                                if 1==chan
                                    plotColor = 'red';
                                else
                                    plotColor = 'green';
                                end                 
                i = 0;

                nSteps = ceil((End - Start)/Step) + 1;

%                 Max_Lr = zeros(sum(cell2mat(cellfun(@length, obj.ROICoordinates, 'uniformoutput', false))), obj.Nchannels); % Assuming the first cell has the same number of channels as the rest
%                 Max_r = zeros(sum(cell2mat(cellfun(@length, obj.ROICoordinates, 'uniformoutput', false))), obj.Nchannels);
%                 Lr_r_Result = zeros(nSteps, sum(cell2mat(cellfun(@length, obj.ROICoordinates, 'uniformoutput', false))), obj.Nchannels);

                Max_Lr = zeros(sum(cell2mat(cellfun(@length, obj.ROICoordinates, 'uniformoutput', false))), 1); % Assuming the first cell has the same number of channels as the rest
                Max_r = zeros(sum(cell2mat(cellfun(@length, obj.ROICoordinates, 'uniformoutput', false))), 1);
                Lr_r_Result = zeros(nSteps, sum(cell2mat(cellfun(@length, obj.ROICoordinates, 'uniformoutput', false))), 1);

                    for roiIter = 1:length(obj.ROICoordinates) % ROI number
                        
                            disp([roiIter length(obj.ROICoordinates)]);
                            
                            q = roiIter;
                            p = 1;

                            roi = obj.ROICoordinates{roiIter};
                            CurrentROI = roi;
                            CurrentROI = [CurrentROI(1,1),  CurrentROI(1,2), max(CurrentROI(:,1)) - min(CurrentROI(:,1)), max(CurrentROI(:,2)) - min(CurrentROI(:,2))];
%                             whichPointsInROI = fliplr(dec2bin(obj.CellData{cellIter}(:,obj.NDataColumns + 1)));
%                             whichPointsInROI = whichPointsInROI(:,roiIter) == '1';
                    x_roi=roi(:,1);
                    y_roi=roi(:,2);
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 

                            dataCropped = obj.CellData{chan}(whichPointsInROI, :);

                            if ~isempty(dataCropped)

                                i=i+1;

                                size_ROI = CurrentROI(3:4);
                                A = polyarea(obj.ROICoordinates{roiIter}(:,1), ...
                                    obj.ROICoordinates{roiIter}(:,2));

                                % Calculate RipleyK for this cell + ROI w/ Channel 1
                                % data

                                assignin('base', 'dataCropped', dataCropped);
                                assignin('base', 'MaxSampledPts', MaxSampledPts);

                                if MaxSampledPts < size(dataCropped)
                                    selectNums = randsample(1:size(dataCropped, 1), MaxSampledPts);
                                    selectVector = false(size(dataCropped, 1), 1);
                                    selectVector(selectNums) = true;
                                else 
                                    selectVector = true(size(dataCropped, 1), 1);
                                end
                                                                    
                                    try
                                    [r, Lr_r] = RipleyKFun(dataCropped(selectVector,5:6), ...
                                                                    A, Start, End, Step, size_ROI);
                                    catch err
                                        disp(err.message);
                                        continue;
                                    end                                                                       
                                    
                                    RipleyKCh1Fig = figure('color', [1 1 1]);
                                    RipleyKCh1Ax = axes('parent', RipleyKCh1Fig);
                                    plot(RipleyKCh1Ax, r, Lr_r, 'color', plotColor, 'linewidth', 2);

                                    % Collect results from these calculations
                                    [MaxLr_r, Index] = max(Lr_r);
                                    Max_Lr(i, 1) = MaxLr_r;     % chan==1
                                    Max_r(i, 1) = r(Index);     % chan==1
                                    Lr_r_Result(:,i, 1) = Lr_r; % chan==1

                                    annotation('textbox', [0.45,0.8,0.22,0.1],...
                                        'String', sprintf('Max L(r) - r: %.3f at Max r : %d', MaxLr_r, Max_r(i)), ...
                                        'FitBoxToText','on');
                                    xlabel(RipleyKCh1Ax, 'r (nm)', 'fontsize', 12);
                                    ylabel(RipleyKCh1Ax, 'L(r) - r', 'fontsize', 12);

                                    print(fullfile(Fun_OutputFolder_name,'RipleyK_Plots',sprintf('Ripley_%dRegion_%d.tif', p, q)), ...
                                        RipleyKCh1Fig, '-dtiff');
                                    close(RipleyKCh1Fig);

                                    Matrix_Result = [r, Lr_r];
                                    SheetName = sprintf('Cell_%dRegion_%d', p, q);
                                    if false %ispc
                                        xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results','RipleyK_Results.xls'), ArrayHeader, SheetName, 'A1');
                                        xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results','RipleyK_Results.xls'), Matrix_Result, SheetName, 'A2');
                                    else
                                        xlsname = fullfile(Fun_OutputFolder_name, 'RipleyK_Results',['RipleyK_Results_ROI_' num2str(roiIter) '.xls']);
                                        save_Excel_or_else(xlsname,ArrayHeader,Matrix_Result); 
                                    end                                        
                            end
                    end

                    % there was loop over channels here..
                    true_chan = chan;
                    chan = 1;
                    
                    Average_Lr_r(:,1) = r;
                    Average_Lr_r(:, 2) = squeeze(mean(Lr_r_Result(:,:,chan), 2));
                    Std_Lr_r(:,2) = std(Lr_r_Result(:,:,chan), 0, 2);

                    Max_r_Ave=[mean(Max_r(:,chan)), std(Max_r(:,chan))];
                    Max_Lr_Ave=[mean(Max_Lr(:,chan)), std(Max_Lr(:,chan))];
                    if false %ispc
                        try                    
                            % Data average on all the regions and cells
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', true_chan)), ArrayHeader, 'Pooled data', 'A1');
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', true_chan)), Average_Lr_r, 'Pooled data', 'A2');

                            % average for max Lr-r and r(max Lr-r)
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', true_chan)), [{'r(max_Lr)'},{'Max_Lr'}], 'Pooled data', 'D3');
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', true_chan)), [{'Mean'},{'Std'}]', 'Pooled data', 'E2');
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', true_chan)), [Max_r_Ave' Max_Lr_Ave'], 'Pooled data', 'E3');

                            % max Lr-r and r(max Lr-r) for each region
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', true_chan)), [{'r(max_Lr)'},{'Max_Lr'}], 'Pooled data', 'E6');
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', true_chan)), [Max_r Max_Lr], 'Pooled data', 'E7');
                        catch err
                            disp(err.message);
                        end
                    else
                        xlsname = fullfile(Fun_OutputFolder_name, 'RipleyK_Results', 'Pooled_Data_Average_Lr_r.xls');
                        save_Excel_or_else(xlsname,ArrayHeader,Average_Lr_r);
                            xlsname = fullfile(Fun_OutputFolder_name, 'RipleyK_Results', 'Pooled_Data_Max_r_Ave_Max_Lr_Ave.xls');
                            captions = [{'r(max_Lr)'},{'Max_Lr'}];
                            datas = [Max_r_Ave' Max_Lr_Ave'];
                            save_Excel_or_else(xlsname,captions,datas);                  
                                xlsname = fullfile(Fun_OutputFolder_name, 'RipleyK_Results', 'Pooled_Data_Max_r_Max_Lr.xls');
                                captions = [{'r(max_Lr)'},{'Max_Lr'}];
                                datas = [Max_r Max_Lr];
                                save_Excel_or_else(xlsname,captions,datas);                            
                            
                    end
                    % use mat format instead
                    Max_r_Ave_Max_Lr_Ave = [Max_r_Ave' Max_Lr_Ave'];
                    Max_r_Max_Lr = [Max_r; Max_Lr]';
                    save(fullfile(Fun_OutputFolder_name, 'RipleyK_Results','Pooled_data.mat'),'Average_Lr_r','Max_r_Ave_Max_Lr_Ave','Max_r_Max_Lr');
    
                    RipleyKMeanFig = figure('color', [1 1 1]);
                    clf(RipleyKMeanFig);
                    RipleyKMeanAx = axes('parent', RipleyKMeanFig, 'nextplot', 'add');
                    plot(RipleyKMeanAx, r, mean(Lr_r_Result(:,:,chan), 2), 'linewidth', 2, 'color', plotColor);
                    plot(RipleyKMeanAx, r, mean(Lr_r_Result(:,:,chan), 2) + Std_Lr_r(:,2), ...
                        'linewidth', 2, 'linestyle', ':', 'color', rgb(52, 152, 219));
                    plot(RipleyKMeanAx, r, mean(Lr_r_Result(:,:,chan), 2) - Std_Lr_r(:,2), ...
                        'linewidth', 2, 'linestyle', ':', 'color', rgb(52, 152, 219));
                    xlabel(RipleyKMeanAx, 'r (nm)', 'fontsize', 12);
                    ylabel(RipleyKMeanAx, 'L(r) - r', 'fontsize', 12);

                    annotation('textbox', [0.45,0.8,0.22,0.1],...
                        'String', sprintf('Max L(r) - r: %.3f at Max r : %d', MaxLr_r, Max_r(i)), ...
                        'FitBoxToText','on');

                    print(fullfile(Fun_OutputFolder_name, 'RipleyK_Plots', 'RipleyK_Average.tif'), ...
                        RipleyKMeanFig, '-dtiff');
                    close(RipleyKMeanFig);
                    % end loop over channels .. ? ..

                valOut = 1;
    end
%-------------------------------------------------------------------------%   
function h = Define_Square_ROIs_Auto(obj,varargin) 
            
            if 1 == nargin
                chan = 1;
            else
                chan = varargin{1};
            end        
            
            Nrois = obj.Square_ROIs_Auto_maxNrois;
            anm = obj.Square_ROIs_Auto_anm;
            nmppix = obj.pixelSizenm;     
            
            % if isnumeric(chan) && intersect(chan,[1 2])
            if strcmp(obj.Square_ROIs_Auto_method,'channel')
                                 
                % to have pixel roughly the size of ROI                 
                size_of_ROI_as_pixel_in_microns = anm/1000;
                ROI_numbers = obj.get_localisation_numbers(chan,size_of_ROI_as_pixel_in_microns);                                
                LNT = obj.Square_ROIs_Auto_LNT(chan);    % low number threshold
                HNT = obj.Square_ROIs_Auto_HNT(chan);    % high number threshold                
                numbers_OK  = LNT < ROI_numbers&ROI_numbers < HNT; 

                qthresh = obj.Square_ROIs_Auto_qthresh(chan);                  
                % to have pixel roughly the size of ROI   
                % first step - just back to widefield                        
%                 ash = obj.get_ash(nmppix,chan);                
%                 f = anm/nmppix;
%                 z = imresize(ash,1/f);          
%                 z = z.*(z > quantile(z(:),qthresh)); 
                z = ROI_numbers;          
                z = z.*(z > quantile(z(:),qthresh));
                %
                obj.ROICoordinates = cell(0);
                for k=1:size(z,1)
                    for m=1:size(z,2)
                        if z(k,m)>0 && numbers_OK(k,m)
                           x = round((m-0.5)*size_of_ROI_as_pixel_in_microns*1000); % back to nanometers to work with table
                           y = round((k-0.5)*size_of_ROI_as_pixel_in_microns*1000);                           
                          d = floor(anm/2)-1;
                          p1 = [x-d y+d];
                          p2 = [x+d y+d];
                          p3 = [x+d y-d];
                          p4 = [x-d y-d];
                          p5 = [x-d y+d];                      
                          roi = [p1;p2;p3;p4;p5];                      
                         obj.ROICoordinates = [obj.ROICoordinates; roi];
                        end
                    end
                end
                
            elseif strcmp(obj.Square_ROIs_Auto_method,'composite')
                
                % to have pixel roughly the size of ROI                 
                size_of_ROI_as_pixel_in_microns = anm/1000;
                ROI_numbers_1 = obj.get_localisation_numbers(1,size_of_ROI_as_pixel_in_microns);                                
                LNT_1 = obj.Square_ROIs_Auto_LNT(1);      % low number threshold
                HNT_1 = obj.Square_ROIs_Auto_HNT(1);    % high number threshold                
                numbers_OK_1  = LNT_1 < ROI_numbers_1&ROI_numbers_1 < HNT_1; 
                
                ROI_numbers_2 = obj.get_localisation_numbers(2,size_of_ROI_as_pixel_in_microns);                                
                LNT_2 = obj.Square_ROIs_Auto_LNT(2);      % low number threshold
                HNT_2 = obj.Square_ROIs_Auto_HNT(2);    % high number threshold               
                numbers_OK_2  = LNT_2 < ROI_numbers_2&ROI_numbers_2 < HNT_2;      
                               
                qthresh_1 = obj.Square_ROIs_Auto_qthresh(1);                  
                qthresh_2 = obj.Square_ROIs_Auto_qthresh(2);
                % to have pixel roughly the size of ROI   
                % first step - just back to widefield                        
%                 ash_1 = obj.get_ash(nmppix,1);
%                 ash_2 = obj.get_ash(nmppix,2);
%                 %
%                 f = anm/nmppix;
%                 z_1 = imresize(ash_1,1/f);          
%                 z_1 = z_1.*(z_1 > quantile(z_1(:),qthresh_1));
%                 %
%                 z_2 = imresize(ash_2,1/f);          
%                 z_2 = z_2.*(z_2 > quantile(z_2(:),qthresh_2));

                z_1 = ROI_numbers_1;          
                z_1 = z_1.*(z_1 > quantile(z_1(:),qthresh_1));
                %
                z_2 = ROI_numbers_2;          
                z_2 = z_2.*(z_2 > quantile(z_2(:),qthresh_2));
                %                
                obj.ROICoordinates = cell(0);
                for k=1:size(z_1,1)
                    for m=1:size(z_1,2)
                        if z_1(k,m)>0 && z_2(k,m)>0 && numbers_OK_1(k,m) && numbers_OK_2(k,m)
                           x = round((m-.5)*size_of_ROI_as_pixel_in_microns*1000); % back to nanometers to work with table
                           y = round((k-.5)*size_of_ROI_as_pixel_in_microns*1000);                           
                          d = floor(anm/2)-1;
                          p1 = [x-d y+d];
                          p2 = [x+d y+d];
                          p3 = [x+d y-d];
                          p4 = [x-d y-d];
                          p5 = [x-d y+d];                      
                          roi = [p1;p2;p3;p4;p5];                      
                         obj.ROICoordinates = [obj.ROICoordinates; roi];
                        end
                    end
                end                                                             
            end

% check rois
bad_rois = [];
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        if strcmp(obj.Square_ROIs_Auto_method,'channel')
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x1=dataCropped(:,5);

                            if numel(x1)<obj.Square_ROIs_Auto_LNT
                                bad_rois = [bad_rois; roiInc];
                            end
                        elseif strcmp(obj.Square_ROIs_Auto_method,'composite')                            
                            chan = 1;
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x1=dataCropped(:,5);
                            chan = 2;
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x2=dataCropped(:,5);                                            

                            if numel(x1)<obj.Square_ROIs_Auto_LNT(1) || numel(x2)<obj.Square_ROIs_Auto_LNT(2)
                                bad_rois = [bad_rois; roiInc];
                            end 
                        end
                  end
%bad_rois
numel(obj.ROICoordinates)
% check rois
            obj.ROICoordinates = obj.ROICoordinates(setdiff(1:numel(obj.ROICoordinates),bad_rois));            
            % if too many ROIs, choose random Nrois among defined
            if (numel(obj.ROICoordinates) > Nrois)
                obj.ROICoordinates = obj.ROICoordinates(randi(numel(obj.ROICoordinates),1,Nrois));
            end             

h = [];            
if obj.Square_ROIs_Auto_verbose
% display
YMAX = obj.SizeY*obj.pixelSizenm;
h = figure('units','normalized','outerposition',[0 0 1 1],'name','ROIs as they go..');
    if strcmp(obj.Square_ROIs_Auto_method,'composite')
                        ax = gca;
                        plot(ax,obj.CellData{1}(:,5),YMAX-obj.CellData{1}(:,6),'r.',obj.CellData{2}(:,5),YMAX-obj.CellData{2}(:,6),'g.');
% markersize = 4;                        
% h = scatter(ax,obj.CellData{1}(:,5),YMAX-obj.CellData{1}(:,6),markersize,'red','filled');
%         alpha = 0.3;
%         set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
% hold on;
% h = scatter(ax,obj.CellData{2}(:,5),YMAX-obj.CellData{2}(:,6),markersize,'green','filled');
%         alpha = 0.3;
%         set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                        
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        chan = 1;
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x1=dataCropped(:,5);
                        y1=dataCropped(:,6);
                        chan = 2;
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x2=dataCropped(:,5);
                        y2=dataCropped(:,6);                        
                        
                        color = [0.2 0.2 0.2+0.79*rand];
                                               
                        plot(ax,x_roi,YMAX-y_roi,'marker','.','color',color,'linestyle','-','linewidth',3);
                        hold(ax,'on');
                        
                        if isempty(x1)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 1]);
                        end 
                        
                        if isempty(x2)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 2]);
                        end                         
                        
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');                  
                   end
                   hold(ax,'off');
                   
    elseif strcmp(obj.Square_ROIs_Auto_method,'channel')
                        ax = gca;
                        if 1==chan
                            clor = 'r.';
                        else
                            clor = 'g.';
                        end
                        plot(ax,obj.CellData{chan}(:,5),YMAX-obj.CellData{chan}(:,6),clor);
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                        
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x1=dataCropped(:,5);
                        
                        color = [0.2 0.2 0.2+0.79*rand];
                                               
                        plot(ax,x_roi,YMAX-y_roi,'marker','.','color',color,'linestyle','-','linewidth',3);
                        hold(ax,'on');
                        
                        if isempty(x1)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 1]);
                        end 
                        
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');                  
                   end
                   hold(ax,'off');  
    end
end % verbose
% display            
                                   
end
%-------------------------------------------------------------------------%
function v = get_localisation_numbers(obj,chan,umppix)                     
                        nmppix = obj.pixelSizenm;
                        x=obj.CellData{chan}(:,5)/1000/umppix;
                        y=obj.CellData{chan}(:,6)/1000/umppix;
                        Lx = round(obj.SizeX*nmppix/1000/umppix);
                        Ly = round(obj.SizeY*nmppix/1000/umppix);
                        x=int64(round(x));
                        y=int64(round(y));                        
                        mask = x>=1 & x<=Lx & y>=1 & y<=Ly;
                        x=x(mask);
                        y=y(mask);
                        v = zeros(Lx,Ly);
                         for k = 1:numel(x)
                            x_k = x(k);
                            y_k = y(k);
                            v(y_k,x_k) = v(y_k,x_k) + 1;
                         end
end
%-------------------------------------------------------------------------%
function [dx2,dy2] = Get_channel2_registration_corrections(obj,~) % translation only
     
    dx2 = 0; dy2 = 0;    
    if 2 ~= numel(obj.CellData), return, end

% FOR TESTING    
%             % cheating - introduce shift to channel 2;
%             dx0 = 45
%             dy0 = -70
%             x = obj.CellData{2}(:,5) + dx0;
%             y = obj.CellData{2}(:,6) + dy0;
%             x(x>obj.SizeX*obj.pixelSizenm) = obj.SizeX*obj.pixelSizenm;
%             x(x<1) = 1;
%             y(y>obj.SizeY*obj.pixelSizenm) = obj.SizeY*obj.pixelSizenm;
%             y(y<1) = 1;                        
%             obj.CellData{2}(:,5) = x;
%             obj.CellData{2}(:,6) = y;
           
            % the alignment is image-based, so first get channel images at pre-defined nmmppix
            nmppix = obj.Align_channels_nmppix; 
            
            ash1 = obj.get_ash(nmppix,1);
            ash2 = obj.get_ash(nmppix,2);            
            %
            % normalize images.. make too bright parts less influential
              s1 = ash1(ash1>0);  
              s2 = ash2(ash2>0);
              t1 = quantile(s1(:),0.995);
              ash1(ash1>t1)=t1;
              t2 = quantile(s2(:),0.995);
              ash2(ash2>t2)=t2;              
            %            
            switch obj.Align_channels_method 
                
                case 'xcorr2'
           
                    snm = 16; % scale, nm
                    s1 = round(max(2,snm/nmppix));
                    s2 = 2*s1;
                    tform = Correct_Simple_Image_Translation_xcorr2(ash2,ash1,s1,s2);

                case 'Matlab_multimodal' 

                        [optimizer, metric] = imregconfig('multimodal');
                     
                        % Tune the properties of the optimizer to get the problem to converge
                        % on a global maxima and to allow for more iterations.
                    %     optimizer.InitialRadius = 0.009;
                    %     optimizer.Epsilon = 1.5e-4;
                    %     optimizer.GrowthFactor = 1.01;
                    %     optimizer.MaximumIterations = 300;
                     
                        tform = imregtform(ash2, ash1, 'translation', optimizer, metric);
                     
            end

% FOR VISUAL ASSESSMENT            
%                     ALYtools_dir = 'C:\Users\alexany\ALYtools';
%                     addpath(ALYtools_dir);
%                     addpath_ALYtools;                          
%                     ash2_reg = imwarp(ash2,tform,'OutputView',imref2d(size(ash1)));
%                     ash = zeros(size(ash1,1),size(ash1,2),3,1,1);
%                     ash(:,:,1,1,1) = ash1;
%                     ash(:,:,2,1,1) = ash2;
%                     ash(:,:,3,1,1) = ash2_reg;
%                     icy_imshow(ash);             
            
             dx2 = - tform.T(3,1)*nmppix;
             dy2 = - tform.T(3,2)*nmppix; % to get it in nanometres
             
             % if shift is bigger than 7 pixels in original image,
             % then the registration is failed
             if norm([dx2 dy2]) > 1000 % OK if not more than micron in either direction
                 disp('registration failed');
                 dx2 = 0;
                 dy2 = 0;
             end
end
%-------------------------------------------------------------------------%
   function ash = get_ash(obj,nmppix,chan,varargin)
              
            X = obj.CellData{chan}(:,5);
            Y = obj.CellData{chan}(:,6);
            
            % if one wants to exclude localisations outside selected ROIs  
            if 4==nargin 
                roi_only_flag = varargin{1};
                if true==roi_only_flag                    
                    whichPointsInROIs = zeros(size(X));
                    for roiInc = 1:length(obj.ROICoordinates)
                        roi = obj.ROICoordinates{roiInc};                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        whichPointsInROI = X>=min(x_roi) & X<=max(x_roi) & Y>=min(y_roi) & Y<=max(y_roi); 
                        whichPointsInROIs = whichPointsInROIs | whichPointsInROI; 
                    end
                    X = X(whichPointsInROIs);
                    Y = Y(whichPointsInROIs);                                        
                end
            end
                         
            SY = ceil(obj.SizeX*obj.pixelSizenm/nmppix);
            SX = ceil(obj.SizeY*obj.pixelSizenm/nmppix);                
            % set diameter to 12 nm
            ash_size = round(max(3,12/nmppix));
            ash = ASH_2d(SX,SY,[Y X]/nmppix,ash_size);
   end
%-------------------------------------------------------------------------%
   function Apply_channel2_registration_corrections(obj,dx2dy2,~)
            X = obj.CellData{2}(:,5) - dx2dy2(1);
            Y = obj.CellData{2}(:,6) - dx2dy2(2);
            nmppix = obj.pixelSizenm;
            mask  = X>=0 & X<=obj.SizeX*nmppix & Y>=0 & Y<=obj.SizeY*nmppix;
            obj.CellData{2} = obj.CellData{2}(mask,:); 
            obj.CellData{2}(:,5) = X(mask);
            obj.CellData{2}(:,6) = Y(mask);            
   end
%-------------------------------------------------------------------------%
    function Analyze_ROIs_DoC(obj,~)
  
        try                      
                   % PLACE COLOCALIZATION RESULTS INTO 1 CHANNEL DIRECTORY
                   dbscanParams = struct;
                       fname = strrep(obj.fileName{1},'.csv','');
                       fname = strrep(fname,'.txt','');
                       DoC_out_dirname = [obj.Outputfolder filesep fname '_channel_' num2str(1) '_ClusDoC_Results' filesep 'DoC'];
                       if ~exist(DoC_out_dirname,'dir')
                           mkdir( fullfile(obj.Outputfolder,[fname '_channel_' num2str(1) '_ClusDoC_Results'],'DoC'));
                           mkdir(DoC_out_dirname,'DBSCAN_Results');
                            mkdir([DoC_out_dirname filesep 'DBSCAN_Results'],'Ch1');
                            mkdir([DoC_out_dirname filesep 'DBSCAN_Results'],'Ch2');
                            mkdir([DoC_out_dirname filesep 'DBSCAN_Results' filesep 'Ch1'],'Cluster_maps');
                            mkdir([DoC_out_dirname filesep 'DBSCAN_Results' filesep 'Ch2'],'Cluster_maps');
                            mkdir([DoC_out_dirname filesep 'DBSCAN_Results' filesep 'Ch1'],'Cluster_density_maps');
                            mkdir([DoC_out_dirname filesep 'DBSCAN_Results' filesep 'Ch2'],'Cluster_density_maps');
                            mkdir(DoC_out_dirname,'DoC_histograms');
                            mkdir(DoC_out_dirname,'DoC_Statistics_and_Plots');
                            mkdir([DoC_out_dirname filesep 'DoC_Statistics_and_Plots'],'Density_and_DoC_maps');
                            mkdir([DoC_out_dirname filesep 'DoC_Statistics_and_Plots'],'Raw_data_maps');
                            mkdir([DoC_out_dirname filesep 'DoC_Statistics_and_Plots'],'Raw_data_maps_with_outliers_removed');
                       end

            dbscanParams(1).Outputfolder = DoC_out_dirname;
            dbscanParams(2).Outputfolder = DoC_out_dirname;
            
            dbscanParams(1).DoCThreshold = obj.DoC_dbscanParams_ch1.DoCThreshold;            
            dbscanParams(1).epsilon = obj.DoC_dbscanParams_ch1.epsilon;
            dbscanParams(1).minPts = obj.DoC_dbscanParams_ch1.minPts;             
            dbscanParams(1).UseLr_rThresh = obj.DoC_dbscanParams_ch1.UseLr_rThresh;
            dbscanParams(1).Lr_rThreshRad = obj.DoC_dbscanParams_ch1.Lr_rThreshRad;
            dbscanParams(1).SmoothingRad = obj.DoC_dbscanParams_ch1.SmoothingRad;
            dbscanParams(1).Cutoff = obj.DoC_dbscanParams_ch1.Cutoff;
            dbscanParams(1).threads = obj.DoC_dbscanParams_ch1.threads;
            dbscanParams(1).DoStats = obj.DoC_dbscanParams_ch1.DoStats;          
            
            dbscanParams(2).DoCThreshold = obj.DoC_dbscanParams_ch2.DoCThreshold;
            dbscanParams(2).epsilon = obj.DoC_dbscanParams_ch2.epsilon;
            dbscanParams(2).minPts = obj.DoC_dbscanParams_ch2.minPts;             
            dbscanParams(2).UseLr_rThresh = obj.DoC_dbscanParams_ch2.UseLr_rThresh;
            dbscanParams(2).Lr_rThreshRad = obj.DoC_dbscanParams_ch2.Lr_rThreshRad;
            dbscanParams(2).SmoothingRad = obj.DoC_dbscanParams_ch2.SmoothingRad;
            dbscanParams(2).Cutoff = obj.DoC_dbscanParams_ch2.Cutoff;
            dbscanParams(2).threads = obj.DoC_dbscanParams_ch2.threads;
            dbscanParams(2).DoStats = obj.DoC_dbscanParams_ch2.DoStats;
                          
            %    original ClusDoC comment:
            % Input parameters:
            % Lr_rad - radius for Lr thresholding check - 20 default
            % Rmax - max distance for DoC Calc (nm) - 500 default
            % Step - step size for DoC Calc (nm) - 10 default
            % ColoThres - threshold for DoC/notDoC - 0.4 default
            % Nb - Num particles with DoC score above threshold to be a 'colocalised' cluster
            % DBSCAN_Radius=20 - epsilon
            % DBSCAN_Nb_Neighbor=3; - minPts ;
            % threads = 2     
            
            Lr_rRad = obj.DoC_Lr_rRad;
            Rmax = obj.DoC_Rmax;
            Step = obj.DoC_Step;
            
            Chan1Color = [.7 0.5 0.3];
            Chan2Color = [.3 0.5 0.7];            
            
            CellData_1 = obj.CellData{1};
            CellData_2 = obj.CellData{2};
            
% for debugging - introduce shift
%             x1 = CellData_1(:,5);
%             y1 = CellData_1(:,6); 
%             y2 = y1 + 70;
%             x2 = x1 + 86;                        
%             mask = x2>0 & x2<=obj.SizeX*obj.pixelSizenm & y2>0 & y2<=obj.SizeY*obj.pixelSizenm;             
%             x2(mask==0) = x1(mask==0);
%             y2(mask==0) = y1(mask==0);             
%             CellData_2(:,5) = x2;
%             CellData_2(:,6) = y2;
%             x1 = CellData_1(:,5);
%             y1 = CellData_1(:,6);
%             x2 = CellData_2(:,5);
%             y2 = CellData_2(:,6);               
%             figure;plot(x1,y1,'r.',x2,y2,'b.');daspect([1 1 1]);            
% for debugging - introduce noise            

             CellData_1(:,12)=1;
             CellData_2(:,12)=2;
             CellData = cat(1,CellData_1,CellData_2);                                        
                
            [DoC_out_CellData, DensityROI] = DoCHandler_YA( ...
                obj.ROICoordinates, ...
                CellData, ... 
                Lr_rRad, ...
                Rmax, ...
                Step, ...
                Chan1Color, ... 
                Chan2Color, ...
                dbscanParams(1).Outputfolder, ...
                obj.NDataColumns);            
                  
            ColoThres = obj.DoC_ColoThres;
            
            ResultTable = ProcessDoCResults_YA( ...
            DoC_out_CellData, ... 
            obj.NDataColumns, ...
            obj.ROICoordinates, ...
            DensityROI, ... 
            DoC_out_dirname, ...
            ColoThres);        

%             % Run DBSCAN on data used for DoC analysis
%             [ClusterTableCh1, ClusterTableCh2, clusterIDOut, handles.ClusterTable] = DBSCANonDoCResults(handles.CellData, handles.ROICoordinates, ...
%                 strcat(handles.Outputfolder, '\Clus-DoC Results'), handles.Chan1Color, handles.Chan2Color, dbscanParams, handles.NDataColumns);

            dbscanParams(1).Outputfolder = [DoC_out_dirname filesep 'DBSCAN_Results' filesep 'Ch1'];
            dbscanParams(2).Outputfolder = [DoC_out_dirname filesep 'DBSCAN_Results' filesep 'Ch2']; % ths won't be reused;
            %
            [ClusterTableCh1, ClusterTableCh2, clusterIDOut, ClusterTable] = obj.DBSCANonDoCResults_YA( ...
                DoC_out_CellData, ...
                [DoC_out_dirname filesep 'DBSCAN_Results'], ...
                Chan1Color, ...
                Chan2Color, ...
                dbscanParams);
            
            DoC_out_CellData = obj.AssignDoCDataToPoints_YA(DoC_out_CellData, clusterIDOut);
            
            NbThresh = obj.DoC_NbThresh; % minimal number of co-localized points            
            EvalStatisticsOnDBSCANandDoCResults_YA(ClusterTableCh1, 1, DoC_out_dirname, NbThresh);
            EvalStatisticsOnDBSCANandDoCResults_YA(ClusterTableCh2, 2, DoC_out_dirname, NbThresh);
    
            % too large
            % obj.ExportDoCDataToCSV(DoC_out_CellData,ClusterTable,DoC_out_dirname); % LOC are localisations 
            
            % not needed, as these ones are saved within "DBSCANonDoCResults_YA" function
            %save([DoC_out_dirname filesep 'ClusterTables.mat'],'ClusterTableCh1','ClusterTableCh2','-v7.3');
                        
        catch mError
            
            disp('ClusDoC processing exited with errors');
            %rethrow(mError);
            disp(mError.message);
            
        end        
    end        
%-------------------------------------------------------------------------%
function DoC_out_CellData = AssignDoCDataToPoints_YA(obj,DoC_in_CellData,clusterIDOut,~)
    
    DoC_out_CellData = DoC_in_CellData;
    
        for roiInc = 1:length(obj.ROICoordinates)

            roi = obj.ROICoordinates{roiInc};                        
            x_roi=roi(:,1);
            y_roi=roi(:,2);
            x=DoC_in_CellData(:,5);
            y=DoC_in_CellData(:,6);            
                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi);
                
            DoC_out_CellData(whichPointsInROI & DoC_out_CellData(:,12) == 1, obj.NDataColumns+3) = clusterIDOut{roiInc, 1, 1};
            DoC_out_CellData(whichPointsInROI & DoC_out_CellData(:,12) == 2, obj.NDataColumns+3) = clusterIDOut{roiInc, 1, 2};
          
            roi_ind = DoC_out_CellData(:,obj.NDataColumns+1);
            roi_ind(whichPointsInROI) = roiInc;
            DoC_out_CellData(:,14) = roi_ind;

        end
end    
%-------------------------------------------------------------------------%
function fname = Save_original_channel2_data_with_XY_registration_corrections(obj,save_dir,dx2dy2,~)

        fname = [];
    
        nmppix = obj.pixelSizenm;

        if strcmp(char(obj.Original_from_file_header{2}{1}(3)),'"x [pix]"') % WindSTORM
            x_ind = 3;
            y_ind = 4;
            multiplier = 1/nmppix; % corrections are in nanometers but WindSTORM format is in pixels
        elseif strcmp(char(obj.Original_from_file_header{2}{1}(2)),'"x [nm]"') % ThunderSTORM
            x_ind = 2;
            y_ind = 3;
            multiplier = 1;            
        elseif strcmp(char(obj.Original_from_file_header{2}{1}(3)),'"x [nm]"') % another ThunderSTORM
            x_ind = 3;
            y_ind = 4;            
            multiplier = 1;
        else % X3
        end

        data = cell2mat(obj.Original_from_file_data{2});

        x = data(:,x_ind) - dx2dy2(1)*multiplier;
        y = data(:,y_ind) - dx2dy2(2)*multiplier;

        mask = x>=0 & x<=obj.SizeX*nmppix & y>=0 & y<=obj.SizeY*nmppix;

        data = data(mask,:);
        data(:,x_ind) = x(mask);
        data(:,y_ind) = y(mask);

        captions = obj.Original_from_file_header{2}{1};

        cap = [];
        for ci=1:numel(captions)
            elem = captions{ci};
            if ci<numel(captions)
                elem = [elem ','];
            end
            cap = [cap elem];
        end
                                    
        fname = obj.fileName{2};
        fname = strrep(fname,'.csv','');
        fname = strrep(fname,'.txt','');
        fname = [fname '_registered.csv'];

        fullfilename = [save_dir filesep fname];

        fid = fopen( fullfilename, 'w' );
        fprintf( fid, '%s\n', cap);
        fclose(fid);
        dlmwrite(fullfilename,data,'-append','precision','%.4f');
end
%-------------------------------------------------------------------------%
function ExportDoCDataToCSV(obj,LOC,ClusterTable,save_dir,~) % LOC are localisations 

        % 1 Index	
        % 2 FirstFrame	
            % 3 NumFrames	
            % 4 NFramesMissing	
        % 5 PostX[nm]	
        % 6 PostY[nm]	
        % 7 Precision[nm]	
        % 8 NPhotons	
        % 9 BkgdVar	
        % 10 Chi^2	
        % 11 PSFWidth[nm]	
        % 12 Channel	
            % 13 ZSlice	
        % 14 ROINum	
            % 15 InOutMask	
        % 16 ClusterID	
        % 17 DoCScore	
        % 18 LrValue	
        % 19 CrossChanDensity	
        % 20 LrAboveThreshold	
        % 21 AllChanDensity

        captions = { ...
        '"Index"', ... % 1 Index	
        '"FirstFrame"', ...% 2 FirstFrame	
        '"NumFrames"', ...    % 3 NumFrames	
        '"NFramesMissing"', ...    % 4 NFramesMissing	
        '"PostX[nm]"', ...% 5 PostX[nm]	
        '"PostY[nm]"', ...% 6 PostY[nm]	
        '"Precision[nm]"', ...% 7 Precision[nm]	
        '"NPhotons"', ...% 8 NPhotons	
        '"BkgdVar"', ...% 9 BkgdVar	
        '"Chi^2"', ...% 10 Chi^2	
        '"PSFWidth[nm]"', ...% 11 PSFWidth[nm]	
        '"Channel"', ...% 12 Channel	
        '"ZSlice"', ...   % 13 ZSlice	
        '"ROINum"', ...% 14 ROINum	
        '"InOutMask"', ...    % 15 InOutMask	
        '"ClusterID"', ...% 16 ClusterID	
        '"DoCScore"', ...% 17 DoCScore	
        '"LrValue"', ...% 18 LrValue	
        '"CrossChanDensity"', ...% 19 CrossChanDensity	
        '"LrAboveThreshold"', ...% 20 LrAboveThreshold	
        '"AllChanDensity"', ...% 21 AllChanDensity
        };    

                fname = 'DoC_Export_Localisations.csv';
                fullfilename = [save_dir filesep fname];

                cap = [];
                for ci=1:numel(captions)
                    elem = captions{ci};
                    if ci<numel(captions)
                        elem = [elem ','];
                    end
                    cap = [cap elem];
                end

                fid = fopen( fullfilename, 'w' );
                fprintf( fid, '%s\n', cap);
                fclose(fid);
                dlmwrite(fullfilename,LOC,'-append','precision','%.5f');
                
% clusters.. 

        captions = { ...
        '"CellNum"', ...
        '"ROINum"', ...
        '"Channel"', ...
        '"ClusterID"', ...
        '"NPoints"', ...
        '"Nb"', ...
        '"MeanDoCScore"', ... 
        '"Area"', ...
        '"Circularity"', ...
        '"TotalAreaDensity"', ...
        '"AvRelativeDensity"', ... 
        '"MeanDensity"', ...
        '"Nb_In"', ...
        '"NInMask"', ...
        '"NOutMask"', ...
        '"Elongation"', ...        
        };

        fname = 'DoC_Export_Clusters.csv';
        fullfilename = [save_dir filesep fname];
        
               cap = [];
                for ci=1:numel(captions)
                    elem = captions{ci};
                    if ci<numel(captions)
                        elem = [elem ','];
                    end
                    cap = [cap elem];
                end

                fid = fopen( fullfilename, 'w' );
                fprintf( fid, '%s\n', cap);
                fclose(fid);
                dlmwrite(fullfilename,ClusterTable,'-append','precision','%.5f');                            
end    


%-------------------------------------------------------------------------%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = Define_Rectangular_ROIs_from_file(obj,fullfilename,chan,super_res_image_pixelSizenm) 

% for rectangular images - [SX SY] are to be SWAPPED re FIJI, ICY etc.    
%     ac.SizeX = 74;   % 125 in Icy 
%     ac.SizeY = 125;  % 74  in Icy  
    
        Nrois = obj.Square_ROIs_Auto_maxNrois;
    
        raw_rois = rectangular_ROIs_from_file(fullfilename);
        
        n_rois = size(raw_rois,1);
        
        obj.ROICoordinates = cell(n_rois,1);
                    
        raw_rois = round(raw_rois*super_res_image_pixelSizenm);
        
        for k=1:n_rois
            x0 = raw_rois(k,1);
            y0 = raw_rois(k,2);
            W = raw_rois(k,3);
            H = raw_rois(k,4);
            %
            p1 = [x0 y0];
            p2 = [x0 y0+H];
            p3 = [x0+W y0+H];
            p4 = [x0+W y0];
            p5 = p1;
            roi_k = [p1;p2;p3;p4;p5]; 
            %
            roi_k(:,2) = round(obj.SizeX*obj.pixelSizenm) - roi_k(:,2);
            %
            obj.ROICoordinates{k} = roi_k;
        end

% check rois
bad_rois = [];
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        if strcmp(obj.Square_ROIs_Auto_method,'channel')
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x1=dataCropped(:,5);

                            if numel(x1)<obj.Square_ROIs_Auto_LNT(chan)
                                bad_rois = [bad_rois; roiInc];
                            end
                        elseif strcmp(obj.Square_ROIs_Auto_method,'composite')                            
                            chan = 1;
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x1=dataCropped(:,5);
                            chan = 2;
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x2=dataCropped(:,5);                                            

                            if numel(x1)<obj.Square_ROIs_Auto_LNT(1) || numel(x2)<obj.Square_ROIs_Auto_LNT(2)
                                bad_rois = [bad_rois; roiInc];
                            end 
                        end
                  end
%bad_rois
numel(obj.ROICoordinates)
% check rois
            obj.ROICoordinates = obj.ROICoordinates(setdiff(1:numel(obj.ROICoordinates),bad_rois));            
            % if too many ROIs, choose random Nrois among defined
            if (numel(obj.ROICoordinates) > Nrois)
                obj.ROICoordinates = obj.ROICoordinates(randi(numel(obj.ROICoordinates),1,Nrois));
            end             

h = [];            
if obj.Square_ROIs_Auto_verbose
% display
YMAX = obj.SizeY*obj.pixelSizenm;
h = figure('units','normalized','outerposition',[0 0 1 1],'name','ROIs as they go..');
    if strcmp(obj.Square_ROIs_Auto_method,'composite')
                        ax = gca;
                        plot(ax,obj.CellData{1}(:,5),YMAX-obj.CellData{1}(:,6),'r.',obj.CellData{2}(:,5),YMAX-obj.CellData{2}(:,6),'g.');
% markersize = 4;                        
% h = scatter(ax,obj.CellData{1}(:,5),YMAX-obj.CellData{1}(:,6),markersize,'red','filled');
%         alpha = 0.3;
%         set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
% hold on;
% h = scatter(ax,obj.CellData{2}(:,5),YMAX-obj.CellData{2}(:,6),markersize,'green','filled');
%         alpha = 0.3;
%         set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                        
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        chan = 1;
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x1=dataCropped(:,5);
                        y1=dataCropped(:,6);
                        chan = 2;
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x2=dataCropped(:,5);
                        y2=dataCropped(:,6);                        
                        
                        color = [0.2 0.2 0.2+0.79*rand];
                                               
                        plot(ax,x_roi,YMAX-y_roi,'marker','.','color',color,'linestyle','-','linewidth',3);
                        hold(ax,'on');
                        
                        if isempty(x1)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 1]);
                        end 
                        
                        if isempty(x2)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 2]);
                        end                         
                        
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');                  
                   end
                   hold(ax,'off');
                   
    elseif strcmp(obj.Square_ROIs_Auto_method,'channel')
                        ax = gca;
                        if 1==chan
                            clor = 'r.';
                        else
                            clor = 'g.';
                        end
                        plot(ax,obj.CellData{chan}(:,5),YMAX-obj.CellData{chan}(:,6),clor);
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                        
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x1=dataCropped(:,5);
                        
                        color = [0.2 0.2 0.2+0.79*rand];
                                               
                        plot(ax,x_roi,YMAX-y_roi,'marker','.','color',color,'linestyle','-','linewidth',3);
                        hold(ax,'on');
                        
                        if isempty(x1)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 1]);
                        end 
                        
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');                  
                   end
                   hold(ax,'off');  
    end
end % verbose
% display            
                                   
end
%-------------------------------------------------------------------------%
function Load_Segmentation(obj,fulpathtosegmetnationomage)
    obj.sgm = imread(fulpathtosegmetnationomage);
end
%-------------------------------------------------------------------------%
function h = Define_Square_ROIs_From_Segmentation(obj,varargin) 
            
            if 1 == nargin
                chan = 1;
            else
                chan = varargin{1};
            end        
            
            Nrois = obj.Square_ROIs_Auto_maxNrois;
            nmppix = obj.pixelSizenm;     
            anm = obj.Square_ROIs_Auto_anm;

            h = [];
            
            try
                objects = obj.sgm;
                if numel(unique(obj.sgm(:))) == 2 % not labelled
                    objects = bwlabel(obj.sgm);
                end
            catch 
                disp('Define_Square_ROIs_From_Segmentation: wrong or empty segmentation image, cannot continue'); return;
            end
            
            [sx,sy] = size(objects);
            if ~(sx==obj.SizeX && sy==obj.SizeY)
                disp('Define_Square_ROIs_From_Segmentation: wrong segmentation image size, cannot continue'); return;
            end
            
            obj.roi_to_object_lut = [];
                        
            if strcmp(obj.Square_ROIs_Auto_method,'channel')
                                 
                % to have pixel roughly the size of ROI                 
                size_of_ROI_as_pixel_in_microns = anm/1000;
                ROI_numbers = obj.get_localisation_numbers(chan,size_of_ROI_as_pixel_in_microns);                                
                LNT = obj.Square_ROIs_Auto_LNT(chan);    % low number threshold
                HNT = obj.Square_ROIs_Auto_HNT(chan);    % high number threshold                
                numbers_OK  = LNT < ROI_numbers&ROI_numbers < HNT; 
                                    
                z = ROI_numbers;          
                % check every square if it is projected on segmentation           
                proj_sgm = zeros(size(z));
                for k=1:size(z,1)
                    for m=1:size(z,2)
                           x = round((m-.5)*size_of_ROI_as_pixel_in_microns*1000/nmppix); % back to original scale
                           y = round((k-.5)*size_of_ROI_as_pixel_in_microns*1000/nmppix);
                           proj_sgm(k,m) = objects(y,x);
                    end
                end                                                
                z = z.*(0~=proj_sgm);
                %
                obj.ROICoordinates = cell(0);
                for k=1:size(z,1)
                    for m=1:size(z,2)
                        if z(k,m)>0 && numbers_OK(k,m)
                           x = round((m-0.5)*size_of_ROI_as_pixel_in_microns*1000); % back to nanometers to work with table
                           y = round((k-0.5)*size_of_ROI_as_pixel_in_microns*1000);                           
                          d = floor(anm/2)-1;
                          p1 = [x-d y+d];
                          p2 = [x+d y+d];
                          p3 = [x+d y-d];
                          p4 = [x-d y-d];
                          p5 = [x-d y+d];                      
                          roi = [p1;p2;p3;p4;p5];                      
                         obj.ROICoordinates = [obj.ROICoordinates; roi];
                         obj.roi_to_object_lut = [obj.roi_to_object_lut; proj_sgm(k,m)];
                        end
                    end
                end
                
            elseif strcmp(obj.Square_ROIs_Auto_method,'composite')
                
                % to have pixel roughly the size of ROI                 
                size_of_ROI_as_pixel_in_microns = anm/1000;
                ROI_numbers_1 = obj.get_localisation_numbers(1,size_of_ROI_as_pixel_in_microns);                                
                LNT_1 = obj.Square_ROIs_Auto_LNT(1);    % low number threshold
                HNT_1 = obj.Square_ROIs_Auto_HNT(1);    % high number threshold                
                numbers_OK_1  = LNT_1 < ROI_numbers_1&ROI_numbers_1 < HNT_1; 
                
                ROI_numbers_2 = obj.get_localisation_numbers(2,size_of_ROI_as_pixel_in_microns);                                
                LNT_2 = obj.Square_ROIs_Auto_LNT(2);    % low number threshold
                HNT_2 = obj.Square_ROIs_Auto_HNT(2);    % high number threshold               
                numbers_OK_2  = LNT_2 < ROI_numbers_2&ROI_numbers_2 < HNT_2;      
                                              
                numbers_OK  = numbers_OK_1 & numbers_OK_2;

                z = ROI_numbers_1; % doesn't matter
                % check every square if it is projectd on segmentation           
                proj_sgm = zeros(size(z));
                for k=1:size(z,1)
                    for m=1:size(z,2)
                           x = round((m-.5)*size_of_ROI_as_pixel_in_microns*1000/nmppix); % back to original scale
                           y = round((k-.5)*size_of_ROI_as_pixel_in_microns*1000/nmppix);
                           proj_sgm(k,m) = objects(y,x);
                    end
                end                                                
                z = z.*(0~=proj_sgm);                           
                %
                obj.ROICoordinates = cell(0);
                for k=1:size(z,1)
                    for m=1:size(z,2)
                        if z(k,m)>0 && numbers_OK(k,m)
                           x = round((m-0.5)*size_of_ROI_as_pixel_in_microns*1000); % back to nanometers to work with table
                           y = round((k-0.5)*size_of_ROI_as_pixel_in_microns*1000);                           
                          d = floor(anm/2)-1;
                          p1 = [x-d y+d];
                          p2 = [x+d y+d];
                          p3 = [x+d y-d];
                          p4 = [x-d y-d];
                          p5 = [x-d y+d];                      
                          roi = [p1;p2;p3;p4;p5];                      
                         obj.ROICoordinates = [obj.ROICoordinates; roi];
                         obj.roi_to_object_lut = [obj.roi_to_object_lut; proj_sgm(k,m)];                         
                        end
                    end
                end
            end

% check rois
bad_rois = [];
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        if strcmp(obj.Square_ROIs_Auto_method,'channel')
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x1=dataCropped(:,5);

                            if numel(x1)<obj.Square_ROIs_Auto_LNT(chan)
                                bad_rois = [bad_rois; roiInc];
                            end
                        elseif strcmp(obj.Square_ROIs_Auto_method,'composite')                            
                            chan = 1;
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x1=dataCropped(:,5);
                            chan = 2;
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x2=dataCropped(:,5);                                            

                            if numel(x1)<obj.Square_ROIs_Auto_LNT(1) || numel(x2)<obj.Square_ROIs_Auto_LNT(2)
                                bad_rois = [bad_rois; roiInc];
                            end 
                        end
                  end
%bad_rois
numel(obj.ROICoordinates)
% check rois
            obj.ROICoordinates = obj.ROICoordinates(setdiff(1:numel(obj.ROICoordinates),bad_rois));
            obj.roi_to_object_lut = obj.roi_to_object_lut(setdiff(1:numel(obj.ROICoordinates),bad_rois));
            % if too many ROIs, choose random Nrois among defined
            if (numel(obj.ROICoordinates) > Nrois)
                mask = randi(numel(obj.ROICoordinates),1,Nrois);
                obj.ROICoordinates = obj.ROICoordinates(mask);
                obj.roi_to_object_lut = obj.roi_to_object_lut(mask);
            end             

h = [];            
if obj.Square_ROIs_Auto_verbose
% display
YMAX = obj.SizeY*obj.pixelSizenm;
h = figure('units','normalized','outerposition',[0 0 1 1],'name','ROIs as they go..');
    if strcmp(obj.Square_ROIs_Auto_method,'composite')
                        ax = gca;
                        plot(ax,obj.CellData{1}(:,5),YMAX-obj.CellData{1}(:,6),'r.',obj.CellData{2}(:,5),YMAX-obj.CellData{2}(:,6),'g.');
% markersize = 4;                        
% h = scatter(ax,obj.CellData{1}(:,5),YMAX-obj.CellData{1}(:,6),markersize,'red','filled');
%         alpha = 0.3;
%         set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
% hold on;
% h = scatter(ax,obj.CellData{2}(:,5),YMAX-obj.CellData{2}(:,6),markersize,'green','filled');
%         alpha = 0.3;
%         set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                        
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        chan = 1;
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x1=dataCropped(:,5);
                        y1=dataCropped(:,6);
                        chan = 2;
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x2=dataCropped(:,5);
                        y2=dataCropped(:,6);                        
                        
                        color = [0.2 0.2 0.2+0.79*rand];
                                               
                        plot(ax,x_roi,YMAX-y_roi,'marker','.','color',color,'linestyle','-','linewidth',3);
                        hold(ax,'on');
                        
                        if isempty(x1)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 1]);
                        end 
                        
                        if isempty(x2)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 2]);
                        end                         
                        
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');                  
                   end
                   hold(ax,'off');
                   
    elseif strcmp(obj.Square_ROIs_Auto_method,'channel')
                        ax = gca;
                        if 1==chan
                            clor = 'r.';
                        else
                            clor = 'g.';
                        end
                        plot(ax,obj.CellData{chan}(:,5),YMAX-obj.CellData{chan}(:,6),clor);
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                        
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x1=dataCropped(:,5);
                        
                        color = [0.2 0.2 0.2+0.79*rand];
                                               
                        plot(ax,x_roi,YMAX-y_roi,'marker','.','color',color,'linestyle','-','linewidth',3);
                        hold(ax,'on');
                        
                        if isempty(x1)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 1]);
                        end 
                        
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');                  
                   end
                   hold(ax,'off');  
    end
end % verbose
% display                                                       
end
%-------------------------------------------------------------------------%
function analysis_controller = save_setups(obj,~)

    s = obj.fileName{1};
    s = strrep(strrep(s,'.csv',''),'.txt','');
    save_full_fname = [obj.Outputfolder filesep 'ClusDoC_service_' s '.mat' ];
    
    analysis_controller.pixelSizenm = obj.pixelSizenm;
    analysis_controller.WindSTORM_Sigmapix = obj.WindSTORM_Sigmapix;
    analysis_controller.Square_ROIs_Auto_anm = obj.Square_ROIs_Auto_anm;
    analysis_controller.Square_ROIs_Auto_qthresh = obj.Square_ROIs_Auto_qthresh;
    analysis_controller.Square_ROIs_Auto_maxNrois = obj.Square_ROIs_Auto_maxNrois;
    analysis_controller.Square_ROIs_Auto_method = obj.Square_ROIs_Auto_method;
    analysis_controller.Square_ROIs_Auto_LNT = obj.Square_ROIs_Auto_LNT;
    analysis_controller.Square_ROIs_Auto_HNT = obj.Square_ROIs_Auto_HNT ;      
    analysis_controller.Align_channels_nmppix = obj.Align_channels_nmppix;
    analysis_controller.Align_channels_method = obj.Align_channels_method ;       
    analysis_controller.pathName = obj.pathName;
    analysis_controller.fileName = obj.fileName ;
    analysis_controller.NDataColumns = obj.NDataColumns;
    analysis_controller.SizeX = obj.SizeX;
    analysis_controller.SizeY = obj.SizeY;
    analysis_controller.Nchannels = obj.Nchannels;
    analysis_controller.ROICoordinates = obj.ROICoordinates;
    analysis_controller.Outputfolder = obj.Outputfolder;   
    analysis_controller.sgm = obj.sgm;
    analysis_controller.sgm_base = obj.sgm_base;
    analysis_controller.roi_to_object_lut = obj.roi_to_object_lut;
    analysis_controller.DoC_Lr_rRad = obj.DoC_Lr_rRad;
    analysis_controller.DoC_Rmax = obj.DoC_Rmax;
    analysis_controller.DoC_Step = obj.DoC_Step;
    analysis_controller.DoC_ColoThres = obj.DoC_ColoThres;
    analysis_controller.DoC_NbThresh = obj.DoC_NbThresh;
    %
    analysis_controller.DoC_dbscanParams_ch1 = obj.DoC_dbscanParams_ch1;
    analysis_controller.DoC_dbscanParams_ch2 = obj.DoC_dbscanParams_ch2;
    analysis_controller.DBSCAN = obj.DBSCAN;
    analysis_controller.RipleyK = obj.RipleyK;
    analysis_controller.Plate = obj.Plate;
    analysis_controller.Plate = obj.Well;
    analysis_controller.Plate = obj.Condition;
        
    save(save_full_fname,'analysis_controller');

end

%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
function h = Define_Square_ROIs_From_Compatible_Polygon(obj,object_index,ROI,varargin) 
            
            if 1 == nargin
                chan = 1;
            else
                chan = varargin{1};
            end        
            
            Nrois = obj.Square_ROIs_Auto_maxNrois;
            nmppix = obj.pixelSizenm;     
            anm = obj.Square_ROIs_Auto_anm;

            h = [];
                        
            obj.roi_to_object_lut = [];
                        
            SX = obj.SizeX*obj.pixelSizenm;
            SY = obj.SizeY*obj.pixelSizenm;            
            Nx = floor(SX/anm);
            Ny = floor(SY/anm);
            
            if strcmp(obj.Square_ROIs_Auto_method,'channel')
                                 
                LNT = obj.Square_ROIs_Auto_LNT(chan);    % low number threshold
                HNT = obj.Square_ROIs_Auto_HNT(chan);    % high number threshold                
                                    
                obj.ROICoordinates = cell(0);
                        
                x=obj.CellData{chan}(:,5);
                y=obj.CellData{chan}(:,6);
                
                for six = 1:Nx
                    for siy = 1:Ny
                        x0 = 1 + anm*(six-1);
                        y0 = 1 + anm*(siy-1);
                        p1 = [x0 y0];
                        p2 = [x0 y0+anm];
                        p3 = [x0+anm y0+anm];
                        p4 = [x0+anm y0];
                        p5 = p1;
                        subROI = [p1;p2;p3;p4;p5];                         
                        in = inpolygon(subROI(:,1),subROI(:,2),ROI(:,1),ROI(:,2));
                        %
                        % look if it satisfies numbers of localizations limitations
                        x_roi=subROI(:,1);
                        y_roi=subROI(:,2);
                        
                        if 5~=sum(in)
                            continue;
                        else
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            n_locs = size(dataCropped,1);
                            if LNT < n_locs&&n_locs < HNT
                             obj.ROICoordinates = [obj.ROICoordinates; subROI];
                             obj.roi_to_object_lut = [obj.roi_to_object_lut; object_index];                            
                            end
                        end 
                    end
                end
            elseif strcmp(obj.Square_ROIs_Auto_method,'composite')

                x1=obj.CellData{1}(:,5);
                y1=obj.CellData{1}(:,6);
                x2=obj.CellData{2}(:,5);
                y2=obj.CellData{2}(:,6);                
                
                LNT_1 = obj.Square_ROIs_Auto_LNT(1);    % low number threshold
                HNT_1 = obj.Square_ROIs_Auto_HNT(1);    % high number threshold                
                
                LNT_2 = obj.Square_ROIs_Auto_LNT(2);    % low number threshold
                HNT_2 = obj.Square_ROIs_Auto_HNT(2);    % high number threshold               

                for six = 1:Nx
                    for siy = 1:Ny
                        x0 = 1 + anm*(six-1);
                        y0 = 1 + anm*(siy-1);
                        p1 = [x0 y0];
                        p2 = [x0 y0+anm];
                        p3 = [x0+anm y0+anm];
                        p4 = [x0+anm y0];
                        p5 = p1;
                        subROI = [p1;p2;p3;p4;p5];                         
                        in = inpolygon(subROI(:,1),subROI(:,2),ROI(:,1),ROI(:,2));
                        % look if it satisfies numbers of localizations limitations
                        x_roi=subROI(:,1);
                        y_roi=subROI(:,2);                        
                        if 5~=sum(in)
                            continue;
                        else
                            whichPointsInROI_1 = x1>=min(x_roi) & x1<=max(x_roi) & y1>=min(y_roi) & y1<=max(y_roi); 
                            dataCropped_1 = obj.CellData{chan}(whichPointsInROI_1,:);
                            n_locs_1 = size(dataCropped_1,1);
                            whichPointsInROI_2 = x2>=min(x_roi) & x2<=max(x_roi) & y2>=min(y_roi) & y2<=max(y_roi); 
                            dataCropped_2 = obj.CellData{chan}(whichPointsInROI_2,:);
                            n_locs_2 = size(dataCropped_2,1);
                            
                            if (LNT_1 < n_locs_1&&n_locs_1 < HNT_1) && (LNT_2 < n_locs_2&&n_locs_2 < HNT_2)
                             obj.ROICoordinates = [obj.ROICoordinates; subROI];
                             obj.roi_to_object_lut = [obj.roi_to_object_lut; object_index];                            
                            end
                        end 
                    end
                end   
            end

% check rois
bad_rois = [];
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        if strcmp(obj.Square_ROIs_Auto_method,'channel')
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x1=dataCropped(:,5);

                            if numel(x1)<obj.Square_ROIs_Auto_LNT(chan)
                                bad_rois = [bad_rois; roiInc];
                            end
                        elseif strcmp(obj.Square_ROIs_Auto_method,'composite')                            
                            chan = 1;
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x1=dataCropped(:,5);
                            chan = 2;
                            x=obj.CellData{chan}(:,5);
                            y=obj.CellData{chan}(:,6);
                                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                            dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                            x2=dataCropped(:,5);                                            

                            if numel(x1)<obj.Square_ROIs_Auto_LNT(1) || numel(x2)<obj.Square_ROIs_Auto_LNT(2)
                                bad_rois = [bad_rois; roiInc];
                            end 
                        end
                  end
%bad_rois
numel(obj.ROICoordinates)
% check rois
            obj.ROICoordinates = obj.ROICoordinates(setdiff(1:numel(obj.ROICoordinates),bad_rois));
            obj.roi_to_object_lut = obj.roi_to_object_lut(setdiff(1:numel(obj.ROICoordinates),bad_rois));
            % if too many ROIs, choose random Nrois among defined
            if (numel(obj.ROICoordinates) > Nrois)
                mask = randi(numel(obj.ROICoordinates),1,Nrois);
                obj.ROICoordinates = obj.ROICoordinates(mask);
                obj.roi_to_object_lut = obj.roi_to_object_lut(mask);
            end             

h = [];            
if obj.Square_ROIs_Auto_verbose
% display
YMAX = obj.SizeY*obj.pixelSizenm;
h = figure('units','normalized','outerposition',[0 0 1 1],'name','ROIs as they go..');
    if strcmp(obj.Square_ROIs_Auto_method,'composite')
                        ax = gca;
                        plot(ax,obj.CellData{1}(:,5),YMAX-obj.CellData{1}(:,6),'r.',obj.CellData{2}(:,5),YMAX-obj.CellData{2}(:,6),'g.');
% markersize = 4;                        
% h = scatter(ax,obj.CellData{1}(:,5),YMAX-obj.CellData{1}(:,6),markersize,'red','filled');
%         alpha = 0.3;
%         set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
% hold on;
% h = scatter(ax,obj.CellData{2}(:,5),YMAX-obj.CellData{2}(:,6),markersize,'green','filled');
%         alpha = 0.3;
%         set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                        
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        chan = 1;
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x1=dataCropped(:,5);
                        y1=dataCropped(:,6);
                        chan = 2;
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x2=dataCropped(:,5);
                        y2=dataCropped(:,6);                        
                        
                        color = [0.2 0.2 0.2+0.79*rand];
                                               
                        plot(ax,x_roi,YMAX-y_roi,'marker','.','color',color,'linestyle','-','linewidth',3);
                        hold(ax,'on');
                        
                        if isempty(x1)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 1]);
                        end 
                        
                        if isempty(x2)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 2]);
                        end                         
                        
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');                  
                   end
                   hold(ax,'off');
                   
    elseif strcmp(obj.Square_ROIs_Auto_method,'channel')
                        ax = gca;
                        if 1==chan
                            clor = 'r.';
                        else
                            clor = 'g.';
                        end
                        plot(ax,obj.CellData{chan}(:,5),YMAX-obj.CellData{chan}(:,6),clor);
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                        
                   for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                        x_roi=roi(:,1);
                        y_roi=roi(:,2);
                        
                        x=obj.CellData{chan}(:,5);
                        y=obj.CellData{chan}(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                        dataCropped = obj.CellData{chan}(whichPointsInROI,:);                       
                        x1=dataCropped(:,5);
                        
                        color = [0.2 0.2 0.2+0.79*rand];
                                               
                        plot(ax,x_roi,YMAX-y_roi,'marker','.','color',color,'linestyle','-','linewidth',3);
                        hold(ax,'on');
                        
                        if isempty(x1)
                            plot(ax,x_roi,YMAX-y_roi,'k:','linewidth',7);
                            hold(ax,'on');
                            %disp([roiInc 1]);
                        end 
                        
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');                  
                   end
                   hold(ax,'off');  
    end
    
end % verbose
% display                                                       
end

%-------------------------------------------------------------------------%
    function Analyze_ROIs_gr(obj,~) 
        try   
               for chan = 1 : obj.Nchannels                   
                   fname = strrep(obj.fileName{chan},'.csv','');
                   fname = strrep(fname,'.txt','');
                   gr_out_dirname = [obj.Outputfolder filesep fname '_channel_' num2str(chan) '_ClusDoC_Results' filesep  'gr'];
                   if ~exist(gr_out_dirname,'dir')
                       mkdir( fullfile(obj.Outputfolder,[fname '_channel_' num2str(chan) '_ClusDoC_Results'],'gr'));
                       mkdir(gr_out_dirname,'gr_Plots');
                       mkdir(gr_out_dirname,'gr_Results');
                   end                                      
                   [~] = obj.grHandler(gr_out_dirname,chan);            
               end                                   
        catch mError            
            disp('g(r) processing exited with errors');
            disp(mError.message);             
        end        
    end
    %-------------------------------------------------------------------------%   
    function valOut = grHandler(obj,Fun_OutputFolder_name,chan,~)
        
        if isempty(obj.ROICoordinates), return, end
        
                ArrayHeader = [{'r'},{'g(r)'}];
                
                                if 1==chan
                                    plotColor = 'red';
                                else
                                    plotColor = 'green';
                                end                 

        Start = obj.RipleyK.Start;
        End = obj.RipleyK.End;
        Step = obj.RipleyK.Step;
        
        nSteps_FFT = obj.RipleyK.End;     
                                                        
                    for roiIter = 1:length(obj.ROICoordinates) % ROI number
                        
                            disp([roiIter length(obj.ROICoordinates)]);
                            
                            roi = obj.ROICoordinates{roiIter};
                            CurrentROI = roi;
                            CurrentROI = [CurrentROI(1,1),  CurrentROI(1,2), max(CurrentROI(:,1)) - min(CurrentROI(:,1)), max(CurrentROI(:,2)) - min(CurrentROI(:,2))];

                            x_roi=roi(:,1);
                            y_roi=roi(:,2);
                                x=obj.CellData{chan}(:,5);
                                y=obj.CellData{chan}(:,6);
                                    whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 

                            dataCropped = obj.CellData{chan}(whichPointsInROI, :);

                            if ~isempty(dataCropped)
                                
                                size_ROI = CurrentROI(3:4);
                                A = polyarea(obj.ROICoordinates{roiIter}(:,1), ...
                                    obj.ROICoordinates{roiIter}(:,2));
                                
                                MaxSampledPts = obj.RipleyK.MaxSampledPts;
                                
                                assignin('base', 'dataCropped', dataCropped);
                                assignin('base', 'MaxSampledPts', MaxSampledPts);
                                if MaxSampledPts < size(dataCropped)
                                    selectNums = randsample(1:size(dataCropped, 1), MaxSampledPts);
                                    selectVector = false(size(dataCropped, 1), 1);
                                    selectVector(selectNums) = true;
                                else 
                                    selectVector = true(size(dataCropped, 1), 1);
                                end
                                                                
                             coord = dataCropped(:,5:6);
                             N_locs = size(coord,1);
                             [r_FFT, gr_FFT] = grFunFFT(coord,CurrentROI,nSteps_FFT); 
                             coord = dataCropped(selectVector,5:6);
                             [r_Ripley, gr_Ripley] = grFun(coord,A,Start,End,Step,size_ROI);                             
                             h = figure;
                             ax=gca;
                             plot(ax,r_Ripley,gr_Ripley,'r.-',r_FFT,gr_FFT,'b.-','markersize',12);
                             grid(ax,'on');
                             xlabel(ax,'distance [nm]','fontsize', 12);
                             ylabel(ax,'g(r)','fontsize', 12);
                             legend(ax,{'Ripley','FFT'});
                             figsavename = [Fun_OutputFolder_name filesep 'gr_Plots' filesep 'gr_Region_' num2str(roiIter) '.fig'];
                             savefig(h,figsavename); 
                             close(h);
                             %
                             xlssavename_Ripley = [Fun_OutputFolder_name filesep 'gr_Results' filesep 'gr_Ripley_Region_' num2str(roiIter) '.xls'];
                             save_Excel_or_else(xlssavename_Ripley,ArrayHeader,[r_Ripley, gr_Ripley]);
                             xlssavename_FFT    = [Fun_OutputFolder_name filesep 'gr_Results' filesep 'gr_FFT_Region_' num2str(roiIter) '.xls'];
                             save_Excel_or_else(xlssavename_FFT,ArrayHeader,[r_FFT, gr_FFT]); 
                             %
                             xlsinfosavename = [Fun_OutputFolder_name filesep 'gr_Results' filesep 'Info_' num2str(roiIter) '.xls'];
                             save_Excel_or_else(xlsinfosavename,[{'N_locs'},{'Area'}],[N_locs, A]);                              
                             %                               
                            end
                    end
                valOut = 1;
        %------------------------------------------------        
        function [r, gr] = grFunFFT(coord,CurrentROI,rmax)
            
            D = obj.Square_ROIs_Auto_anm;
            
            X = coord(:,1);
            Y = coord(:,2);
            X0 = CurrentROI(1);
            Y0 = CurrentROI(2) - D;

            u = zeros(D);
            %
                for k=1:length(X)
                        cur_x = round(X(k))-round(X0);
                        cur_y = round(Y(k))-round(Y0);
                        if cur_x>0 && cur_x<=D && cur_y>0 && cur_y<=D
                            u(cur_x,cur_y) = u(cur_x,cur_y)+1;
                        end
                end 
                
            [~, r, gr, ~, ~] = get_autocorr(u , ones(size(u)), rmax, false);
    
            r = r(2:(rmax+1))';
            gr = gr(2:(rmax+1))';                              
        end                
    end
   
%-------------------------------------------------------------------------%
    function gr_perObject_from_Segmentation_Handler(obj,dir_name,rmax,chan,~)
                                if 1==chan
                                    plotColor = 'red';
                                else
                                    plotColor = 'green';
                                end         
        s = regionprops(obj.sgm,'BoundingBox');
        
        bb = round(s(1).BoundingBox);          
        x0 = bb(1);
        y0 = bb(2);
        W = bb(3);
        H = bb(4);
        
        u1 = obj.sgm==1;
      
        bb_SR   = round(s(1).BoundingBox*obj.pixelSizenm);
        x0_SR   = bb_SR(1);
        y0_SR   = bb_SR(2);
        W_SR    = bb_SR(3);
        H_SR    = bb_SR(4);

        p1 = [x0_SR         y0_SR];
        p2 = [x0_SR         y0_SR+H_SR];
        p3 = [x0_SR+W_SR    y0_SR+H_SR];
        p4 = [x0_SR+W_SR    y0_SR];
        p5 = p1;
        roi = [p1;p2;p3;p4;p5];  
        
        x=obj.CellData{chan}(:,5);
        y=obj.CellData{chan}(:,6);                                
        %          
        x_roi=roi(:,1);
        y_roi=roi(:,2);            
        whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi);

        dataCropped = round(obj.CellData{chan}(whichPointsInROI,:));        
        
%         xmax = min(x0+W,obj.SizeX);
%         ymax = min(y0+H,obj.SizeY); % ?? worked because always square images?
        xmax = min(x0+W,obj.SizeY);
        ymax = min(y0+H,obj.SizeX);
        u0 = uint8(imresize(u1(y0:ymax,x0:xmax),obj.pixelSizenm));
        
        u = zeros(size(u0),'single');
        x_SR = dataCropped(:,5);
        y_SR = dataCropped(:,6);        
        for m=1:size(dataCropped,1)
            y_m = x_SR(m) - x0_SR;
            x_m = y_SR(m) - y0_SR;
            try            
                if 0~= u0(x_m,y_m)
                    u(x_m,y_m) = u(x_m,y_m) + 1;
                end
            catch
%                 disp('error')
%                 [x_m y_m]
            end            
        end
        
        % save Object information
        save_Excel_or_else([dir_name filesep 'Object_info.xls'], ...
            [{'N_locs'},{'Area'}],[sum(u,'All'), sum(u0,'All')]);        
        % save binary 1 nm map
        imwrite(uint8(u0&~u),[dir_name filesep 'Object_nm_map.tif']);
        
        tic        
        [~, r, gr, ~, ~] = get_autocorr(single(u) , single(u0), rmax, false);
        r = r(2:(rmax+1))';
        gr = gr(2:(rmax+1))';
        toc/60 

        xlsname = [dir_name filesep 'perOjbect_gr.xls'];
        save_Excel_or_else(xlsname,[{'r'},{'g(r)'}],[r, gr]);
        figname = [dir_name filesep 'perOjbect_gr.fig'];
        h = figure;
        ax=gca;
        plot(ax,r,gr,'color', plotColor, 'linewidth', 2, 'marker','o','markersize',4);
        xlabel(ax,'distance [nm]','fontsize',14);
        ylabel(ax,'AU','fontsize',14);
        grid(ax,'on');
        savefig(h,figname);
        close(h);

if ~obj.GR_ASSISTED_CLUSTERING.PERFORM, return, end

r_start = obj.GR_ASSISTED_CLUSTERING.r_start;
dr = obj.GR_ASSISTED_CLUSTERING.dr;
comp_mode = obj.GR_ASSISTED_CLUSTERING.MODE;
%
fitres = zeros(length(dr),3); % s1,s2,fval
% needed to feed fitting proc
N_locs = sum(u,'All');
Area = sum(u0,'All');
% calculating interval of indices - to find best estimeate for 2 sigmas
r_step = r(2)-r(1);
DN = round(dr/r_step);
for k=1:DN    
    INTERVAL = r<(r_start + r_step*k);
    [g_fit, A, L, n, N, F, fval] = fit_gr_14092023(r(INTERVAL),gr(INTERVAL),N_locs,Area,comp_mode,[]);
    fitres(k,1) = L(1);
    fitres(k,2) = L(2);    
    fitres(k,3) = fval/length(r(INTERVAL));
end
optimal_fit_ind = find(fitres(:,3)==min(fitres(:,3)));
s1 = round(fitres(optimal_fit_ind,1));
s2 = round(fitres(optimal_fit_ind,2));
%
out = obj.gr_assisted_quasi_cluster_analysis(u>0, u0, s1, s2, dir_name);
% g(r)-assisted "true clustering" - end
        
    end
    
%-------------------------------------------------------------------------%
    function Analyze_gr_perObject_from_Segmentation(obj,rmax,~)
  
         try   
               for chan = 1 : obj.Nchannels                   
                   fname = strrep(obj.fileName{chan},'.csv','');
                   fname = strrep(fname,'.txt','');
                   gr_out_dirname = [obj.Outputfolder filesep fname '_channel_' num2str(chan) '_ClusDoC_Results' filesep  'gr'];
                   if ~exist(gr_out_dirname,'dir')
                       mkdir( fullfile(obj.Outputfolder,[fname '_channel_' num2str(chan) '_ClusDoC_Results'],'gr'));
                   end                                      
                   obj.gr_perObject_from_Segmentation_Handler(gr_out_dirname,rmax,chan); 
               end                                   
        catch mError            
            disp('g(r) processing exited with errors');
            disp(mError.message);             
        end        
    end   
%-------------------------------------------------------------------------%    
    end % methods
            
end
