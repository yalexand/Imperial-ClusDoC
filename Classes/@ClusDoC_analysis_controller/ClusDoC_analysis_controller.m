
classdef ClusDoC_analysis_controller < handle 
    
    % 2020 Imperial College London.
    % this module is an adaptation of ClusDoC software
   
    properties(Constant)
        data_settings_filename = 'ClusDoC_settings.xml';
    end
    
    properties(SetObservable = true)
	
    pixelSizenm = 1;
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
        'threads',2, ...
        'DoStats',true, ...
        'Outputfolder',['C:' filesep]);
 
        Square_ROIs_Auto_anm = 1400;
        Square_ROIs_Auto_qthresh = [.5 .5];
        Square_ROIs_Auto_maxNrois = 20; 
        Square_ROIs_Auto_method = 'channel';         
        
        Align_channels_nmppix = 4;
        Align_channels_method = 'Matlab_multimodal';
    
    % Initialize structure to pass values between GUI components
    CellData = {2,1}; % better it rename this as FOVData
    
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
                    
                    obj.SizeX = settings.SizeX;
                    obj.SizeY = settings.SizeY;                           
        
                    obj.Align_channels_nmppix = settings.Align_channels_nmppix;
                    obj.Align_channels_method = settings.Align_channels_method;                    
                                       
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
                    
                    settings.SizeX = obj.SizeX;
                    settings.SizeY = obj.SizeY; 
                    
                    settings.Align_channels_nmppix = obj.Align_channels_nmppix;
                    settings.Align_channels_method = obj.Align_channels_method;                      
                    
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
      
            if goodZENFile || contains(obj.fileName{chan},'.csv')
                importData = ImportFile(fullfile(obj.pathName{chan},obj.fileName{chan}),obj);
                                
                obj.CellData{chan} = [importData zeros(size(importData, 1), 8)];

%                 obj.CellData{k}(any(isnan(obj.CellData{k}), 2), :) = []; % protection against incomplete line writing in ZEN export
                                                                                 % This breaks import for ThunderSTORMConcatenator output                                                                                  % Commenting this out to allow that format to work.
                obj.NDataColumns = size(importData, 2);
                obj.CellData{chan}(:,obj.NDataColumns + 2) = 1; % All data is in mask until set otherwise
                                             
%                 obj.CellData{chan}(any(obj.CellData{chan}(:, 5) > obj.SizeX), : )= [];
%                 obj.CellData{chan}(any(obj.CellData{chan}(:, 6) > obj.SizeY), : )= [];                
%                 obj.CellData{chan}(any(obj.CellData{chan}(:, 5:6) < 0), : )= [];
                
                obj.Nchannels = max(chan,obj.Nchannels); %min([numel(unique(obj.CellData{chan}(:,12))), 2]); % cap import to 2 channels ever
                
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
                                    DBSCANHandler_YA(dataCropped(dataCropped(:,12) == 1, 5:6), obj.DBSCAN, c, roiInc, ...
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
                    ExportDBSCANDataToExcelFiles(cellROIPair, Result, obj.DBSCAN.Outputfolder, 1);
                    catch err
                        disp('error in function ExportDBSCANDataToExcelFiles');
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
                                appendTable = nan(length(ClusterCh), 15);
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
            rethrow(mError);
             
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

                                    chan = 1;
                                    
                                    plotColor = 'red';

                                    try
                                    [r, Lr_r] = RipleyKFun(dataCropped(selectVector & (dataCropped(:,12) == chan),5:6), ...
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
                                    Max_Lr(i, chan) = MaxLr_r;
                                    Max_r(i, chan) = r(Index);
                                    Lr_r_Result(:,i, chan) = Lr_r;

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
                    Average_Lr_r(:,1) = r;
                    Average_Lr_r(:,1) = r;
                    Average_Lr_r(:, 2) = squeeze(mean(Lr_r_Result(:,:,chan), 2));
                    Std_Lr_r(:,2) = std(Lr_r_Result(:,:,chan), 0, 2);

                    Max_r_Ave=[mean(Max_r(:,chan)), std(Max_r(:,chan))];
                    Max_Lr_Ave=[mean(Max_Lr(:,chan)), std(Max_Lr(:,chan))];
                    if false %ispc
                        try                    
                            % Data average on all the regions and cells
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', chan)), ArrayHeader, 'Pooled data', 'A1');
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', chan)), Average_Lr_r, 'Pooled data', 'A2');

                            % average for max Lr-r and r(max Lr-r)
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', chan)), [{'r(max_Lr)'},{'Max_Lr'}], 'Pooled data', 'D3');
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', chan)), [{'Mean'},{'Std'}]', 'Pooled data', 'E2');
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', chan)), [Max_r_Ave' Max_Lr_Ave'], 'Pooled data', 'E3');

                            % max Lr-r and r(max Lr-r) for each region
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', chan)), [{'r(max_Lr)'},{'Max_Lr'}], 'Pooled data', 'E6');
                            xlswrite(fullfile(Fun_OutputFolder_name, 'RipleyK_Results', sprintf('Ch%dPooled.xls', chan)), [Max_r Max_Lr], 'Pooled data', 'E7');
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
function Define_Square_ROIs_Auto(obj,~,varargin) 
            
            if 1 == nargin
                chan = 1;
            else
                chan = varargin(1);
            end        
            
            Nrois = obj.Square_ROIs_Auto_maxNrois;
            anm = obj.Square_ROIs_Auto_anm;
            nmppix = obj.pixelSizenm;     
            
            % if isnumeric(chan) && intersect(chan,[1 2])
            if strcmp(obj.Square_ROIs_Auto_method,'channel')
                 
                qthresh = obj.Square_ROIs_Auto_qthresh(chan);                  

                % to have pixel roughly the size of ROI   
                % first step - just back to widefield                        
                ash = obj.get_ash(nmppix,chan);
                %
                f = anm/nmppix;
                z = imresize(ash,1/f);          
                z = z.*(z > quantile(z(:),qthresh));
                %
                obj.ROICoordinates = cell(0);
                for k=1:size(z,1)
                    for m=1:size(z,2)
                        if z(k,m)>0
                           % disp([k-0.5,m-0.5]);
                           x = round((m-.5)*f*nmppix);
                           y = round((k-.5)*f*nmppix);
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
                
                qthresh_1 = obj.Square_ROIs_Auto_qthresh(1);                  
                qthresh_2 = obj.Square_ROIs_Auto_qthresh(2);
                % to have pixel roughly the size of ROI   
                % first step - just back to widefield                        
                ash_1 = obj.get_ash(nmppix,1);
                ash_2 = obj.get_ash(nmppix,2);
                %
                f = anm/nmppix;
                z_1 = imresize(ash_1,1/f);          
                z_1 = z_1.*(z_1 > quantile(z_1(:),qthresh_1));
                %
                z_2 = imresize(ash_2,1/f);          
                z_2 = z_2.*(z_2 > quantile(z_2(:),qthresh_2));
                %                
                obj.ROICoordinates = cell(0);
                for k=1:size(z_1,1)
                    for m=1:size(z_1,2)
                        if z_1(k,m)>0 && z_2(k,m)>0
                           % disp([k-0.5,m-0.5]);
                           x = round((m-.5)*f*nmppix);
                           y = round((k-.5)*f*nmppix);
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
            % if too many ROIs, choose random Nrois among defined
            if (numel(obj.ROICoordinates) > Nrois)
                obj.ROICoordinates = obj.ROICoordinates(randi(numel(obj.ROICoordinates),1,Nrois));
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
            t1 = quantile(ash1(:),0.995);
            ash1(ash1>t1)=t1;
            t2 = quantile(ash2(:),0.995);
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
             dy2 = - tform.T(3,2)*nmppix;
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
                           mkdir(DoC_out_dirname,'DBSCAN Results');
                            mkdir([DoC_out_dirname filesep 'DBSCAN Results'],'Ch1');
                            mkdir([DoC_out_dirname filesep 'DBSCAN Results'],'Ch2');
                            mkdir([DoC_out_dirname filesep 'DBSCAN Results' filesep 'Ch1'],'Cluster_maps');
                            mkdir([DoC_out_dirname filesep 'DBSCAN Results' filesep 'Ch2'],'Cluster_maps');
                            mkdir([DoC_out_dirname filesep 'DBSCAN Results' filesep 'Ch1'],'Cluster_density_maps');
                            mkdir([DoC_out_dirname filesep 'DBSCAN Results' filesep 'Ch2'],'Cluster_density_maps');
                            mkdir(DoC_out_dirname,'DoC histograms');
                            mkdir(DoC_out_dirname,'DoC Statistics and Plots');
                            mkdir([DoC_out_dirname filesep 'DoC Statistics and Plots'],'Density and DoC maps');
                            mkdir([DoC_out_dirname filesep 'DoC Statistics and Plots'],'Raw data maps');
                            mkdir([DoC_out_dirname filesep 'DoC Statistics and Plots'],'Raw data maps with outliers removed');
                       end

% to setups - all except dirname                       
            dbscanParams(1).Outputfolder = DoC_out_dirname;
            dbscanParams(2).Outputfolder = DoC_out_dirname;            
            dbscanParams(1).DoCThreshold = .4;
            dbscanParams(2).DoCThreshold = .4;
            
            dbscanParams(1).epsilon = 20;
            dbscanParams(1).minPts = 3;             
            dbscanParams(1).UseLr_rThresh = true; 
            dbscanParams(1).Lr_rThreshRad = 20;            
            dbscanParams(1).SmoothingRad = 15;
            dbscanParams(1).Cutoff = 10;            
            dbscanParams(1).threads = 12;
            dbscanParams(1).DoStats = true;            

            dbscanParams(2).epsilon = 20;
            dbscanParams(2).minPts = 3;             
            dbscanParams(2).UseLr_rThresh = true; 
            dbscanParams(2).Lr_rThreshRad = 20;            
            dbscanParams(2).SmoothingRad = 15;
            dbscanParams(2).Cutoff = 10;            
            dbscanParams(2).threads = 12;
            dbscanParams(2).DoStats = true;              
                          
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
            
% to setups  
            Lr_rRad = 20;
            Rmax = 500;
            Step = 10;
            
            Chan1Color = [.7 0.5 0.3];
            Chan2Color = [.3 0.5 0.7];            
            
            CellData_1 = obj.CellData{1};
            CellData_2 = obj.CellData{2};
            
% for debugging - introduce shift
            x1 = CellData_1(:,5);
            y1 = CellData_1(:,6); 
            y2 = y1 + 70;
            x2 = x1 + 86;                        
            mask = x2>0 & x2<=obj.SizeX*obj.pixelSizenm & y2>0 & y2<=obj.SizeY*obj.pixelSizenm;             
            x2(mask==0) = x1(mask==0);
            y2(mask==0) = y1(mask==0);             
            CellData_2(:,5) = x2;
            CellData_2(:,6) = y2;
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
            
% to setups            
            ColoThres = 0.4;
            
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

            dbscanParams(1).Outputfolder = [DoC_out_dirname filesep 'DBSCAN Results' filesep 'Ch1'];
            dbscanParams(2).Outputfolder = [DoC_out_dirname filesep 'DBSCAN Results' filesep 'Ch2']; % ths won't be reused;
            %
            [ClusterTableCh1, ClusterTableCh2, clusterIDOut, handles.ClusterTable] = DBSCANonDoCResults_YA( ...
                DoC_out_CellData, ...
                obj.ROICoordinates, ...
                [DoC_out_dirname filesep 'DBSCAN Results'], ...
                Chan1Color, ...
                Chan2Color, ...
                dbscanParams, ...
                obj.NDataColumns);
                       
            NbThresh = 10;
            EvalStatisticsOnDBSCANandDoCResults_YA(ClusterTableCh1, 1, DoC_out_dirname, NbThresh);
            EvalStatisticsOnDBSCANandDoCResults_YA(ClusterTableCh2, 2, DoC_out_dirname, NbThresh);
            
            save([DoC_out_dirname filesep 'ClusterTables.mat'],'ClusterTableCh1','ClusterTableCh2','-v7.3');
                        
        catch mError
            
            disp('ClusDoC processing exited with errors');
            rethrow(mError);
             
        end        
    end   
             
    end % methods
            
end
