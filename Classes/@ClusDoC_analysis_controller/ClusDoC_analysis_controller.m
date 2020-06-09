
classdef ClusDoC_analysis_controller < handle 
    
    % 2020 Imperial College London.
    % this module is an adaptation of ClusDoC software
   
    properties(Constant)
        data_settings_filename = 'ClusDoC_settings.xml';
    end
    
    properties(SetObservable = true)
	
	src_dir = [];
	dst_dir = [];  
    
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
        Square_ROIs_Auto_qthresh = 20;
        Square_ROIs_Auto_maxNrois = .95;  
    
    % Initialize structure to pass values between GUI components
    CellData = {}; % better it rename this as FOVData
    ROIData = {};
    ROIPos = [];
    CurrentCellData = 1;
    CurrentROIData = [];
    
    pathName = pwd;
    fileName = [];

    % Default ROI settings
    ROISize = 4000; % Length of ROI, in nm   
    
    NDataColumns = [];
    MaxSize = [];
    ROIMultiplier = [];
    Nchannels = []; % ?
    
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
        function Load_Data(obj,fileName,pathName,~)
            
            obj.fileName = fileName;
            obj.pathName = pathName;
        
            goodZENFile = checkZenFile(fullfile(obj.pathName,obj.fileName));
      
            if goodZENFile || contains(obj.fileName,'.csv')
                importData = Import1File(fullfile(obj.pathName,obj.fileName),obj);
                                
                obj.CellData = [importData.Data zeros(size(importData.Data, 1), 8)];
                
%                 obj.CellData{k}(:,5:6) = obj.CellData{k}(:,5:6)*importData.Footer{2}(3)/importData.Footer{2}(1);

%                 obj.CellData{k}(any(isnan(obj.CellData{k}), 2), :) = []; % protection against incomplete line writing in ZEN export
                                                                                 % This breaks import for ThunderSTORMConcatenator output                                                                                  % Commenting this out to allow that format to work.
                obj.NDataColumns = size(importData.Data, 2);
                obj.CellData(:,obj.NDataColumns + 2) = 1; % All data is in mask until set otherwise
                obj.ROIMultiplier = importData.Footer{2}(1); % Conversion from coordinates.txt positions to nm
                             
                % Max size calc has some issues around certain imported ZEN
                % files
                obj.MaxSize = importData.Footer{2}(5)*10*importData.Footer{2}(1)/importData.Footer{2}(3); % FOV size, in nm
                
                if obj.MaxSize == 256
                    obj.MaxSize = obj.MaxSize*100;
                end
                
                % Clear out any points outside of bounds [0 MaxSize];
                obj.CellData(any(obj.CellData(:, 5:6) > obj.MaxSize), : )= [];
                obj.CellData(any(obj.CellData(:, 5:6) < 0), : )= [];
                
                obj.Nchannels = min([numel(unique(obj.CellData(:,12))), 2]); % cap import to 2 channels ever
                
            else
                
                fprintf(1, 'File not in accepted coordinate table format.\nSkipping %s\n', fullfile(obj.pathName,obj.fileName));

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
               chan = 1;
               
               fname = strrep(obj.fileName,'.csv','');
               fname = strrep(fname,'.txt','');
               DBSCAN_out_dirname = [obj.Outputfolder filesep fname '_ClusDoC_Results' filesep 'DBSCAN'];
               if ~exist(DBSCAN_out_dirname,'dir')
                   mkdir( fullfile(obj.Outputfolder,[fname '_ClusDoC_Results'],'DBSCAN'));
                   mkdir(DBSCAN_out_dirname,'Cluster_maps');
                   mkdir(DBSCAN_out_dirname,'Cluster_density_maps');
               end
                          
          obj.DBSCAN.Outputfolder = DBSCAN_out_dirname;  
               
        cellROIPair = [];                  
               
        Result = cell(length(obj.ROICoordinates),1);
        ClusterSmoothTable = cell(length(obj.ROICoordinates),1);
                
                    if verbose
                        figure;  
                        ax = gca;
                        plot(ax,obj.CellData(:,5),obj.CellData(:,6),'b.');
                        daspect(ax,[1 1 1]); 
                        grid(ax,'on');
                        hold(ax,'on');
                    end
                    
                    for roiInc = 1:length(obj.ROICoordinates)

                        roi = obj.ROICoordinates{roiInc};
                        
                    x_roi=roi(:,1);
                    y_roi=roi(:,2);
                        x=obj.CellData(:,5);
                        y=obj.CellData(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                    dataCropped = obj.CellData(whichPointsInROI,:);

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
                            [~, ClusterSmoothTable{roiInc}, ~, classOut, ~, ~, ~, Result{roiInc,c}] = ...
                                DBSCANHandler_YA(dataCropped(dataCropped(:,12) == chan, 5:6), obj.DBSCAN, c, roiInc, ...
                                true, true, clusterColor, dataCropped(dataCropped(:,12) == chan, obj.NDataColumns + 2));

                            % disp(Result);
                            
                            obj.CellData(whichPointsInROI & (obj.CellData(:,12) == chan), obj.NDataColumns + 3) = classOut;
                             
                            cellROIPair = [cellROIPair; c, roiInc, roi(1,1), roi(1,2), polyarea(roi(:,1), roi(:,2))];
                            
                            obj.ClusterTable = AppendToClusterTable(obj.ClusterTable, chan, c, roiInc, ClusterSmoothTable{roiInc, c}, classOut);

                        else
                            % Have chosen an empty region as ROI                            
                            fprintf(1, 'Cell %d - ROI %d is empty.  Skipping.\n', 1, roiInc);                           
                            ClusterSmoothTable{roiInc} = [];
                            classOut = [];
                            Result{roiInc,c} = [];
                        end                        
                    end % ROI
                    
                    if verbose, hold(ax,'off'), end

                if ~all(cellfun(@isempty, Result))
                    ExportDBSCANDataToExcelFiles(cellROIPair, Result, obj.DBSCAN.Outputfolder, chan);
                else
                    fprintf(1, 'All cells and ROIs empty.  Skipping export.\n');
                end
                
               save(fullfile(obj.DBSCAN.Outputfolder,'DBSCAN_Cluster_Result.mat'),'ClusterSmoothTable','Result','-v7.3');                
                                            
        catch mErr
            
            assignin('base', 'cellROIPair', cellROIPair);
            assignin('base', 'Result', Result);
            assignin('base', 'outputFolder', strcat(obj.Outputfolder, '\DBSCAN Results'));
            assignin('base', 'chan', chan);
            
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

               fname = strrep(obj.fileName,'.csv','');
               fname = strrep(fname,'.txt','');
               RipleyK_out_dirname = [obj.Outputfolder filesep fname '_ClusDoC_Results' filesep 'RipleyK'];
               if ~exist(RipleyK_out_dirname,'dir')
                   mkdir( fullfile(obj.Outputfolder,[fname '_ClusDoC_Results'],'RipleyK'));
                   mkdir(RipleyK_out_dirname,'RipleyK_Plots');
                   mkdir(RipleyK_out_dirname,'RipleyK_Results');
               end            
                       
            [~] = obj.RipleyKHandler(RipleyK_out_dirname);
            
        catch mError
            
            disp('Ripley K processing exited with errors');
            rethrow(mError);
             
        end        
    end
    %-------------------------------------------------------------------------%   
    function valOut = RipleyKHandler(obj,Fun_OutputFolder_name,~)
        
        if isempty(obj.ROICoordinates), return, end
        
        Start = obj.RipleyK.Start;
        End = obj.RipleyK.End;
        Step = obj.RipleyK.Step;
        MaxSampledPts = obj.RipleyK.MaxSampledPts;  

            % Handler function for RipleyK calculations
                ArrayHeader = [{'r'},{'L(r)-r'}];

                i = 0;

                nSteps = ceil((End - Start)/Step) + 1;

                Max_Lr = zeros(sum(cell2mat(cellfun(@length, obj.ROICoordinates, 'uniformoutput', false))), obj.Nchannels); % Assuming the first cell has the same number of channels as the rest
                Max_r = zeros(sum(cell2mat(cellfun(@length, obj.ROICoordinates, 'uniformoutput', false))), obj.Nchannels);
                Lr_r_Result = zeros(nSteps, sum(cell2mat(cellfun(@length, obj.ROICoordinates, 'uniformoutput', false))), obj.Nchannels);

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
                        x=obj.CellData(:,5);
                        y=obj.CellData(:,6);
                            whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 

                            dataCropped = obj.CellData(whichPointsInROI, :);

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

                                    [r, Lr_r] = RipleyKFun(dataCropped(selectVector & (dataCropped(:,12) == chan),5:6), ...
                                        A, Start, End, Step, size_ROI);

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
        % Load ROI coordinates from coordinates.txt file (if existing)
        function [roiCoordinates, loadOK] = loadCoordinatesFile(obj,fName,~)

            scaleFactor = obj.ROIMultiplier;
            
            % Optional comment block at top, which may contain line specifying
            % the ROI size.
            % Comments have first character #
            % ROI size specified by # ROISize:\t%f in nanometers
            % If not specified, assume is default value
            % Assuming that all ROIs specified in coordinates.txt file are
            % squares

            fID = fopen(fName, 'r');
            lineNow = 0;
            isEnd = false;
            isData = false;
            while ~isEnd | ~isData
                lineString = fgetl(fID);
                if lineString(1) ~= '#'
                    testLine = lineString;
        %             disp('end of header');
                    isData = true;
                    break;
                elseif isempty(lineString)
                    isEnd = true;
        %             disp('end of file');
                else
                    lineNow = lineNow + 1;
                    % Check if ROI size specified
                    if ~isempty(strfind(lower(lineString), 'roisize'))
                        obj.ROISize = str2double(lineString(regexp(lineString, '\d+'):end));
                    end
                end
            end

            % See if file is ZEN export format of identical rectangles, or is
            % polygons w/ xy coordinates
            % If only 4 columns, then ZEN rectangles
            % Any more than 4 columns and format has to be polygons in
            % x1\ty1\tx2\ty2\tx3\ty3... format

            nTabs = numel(strfind(testLine, sprintf('\t')));
            fseek(fID, 0, -1);
            for skipLines = 1:lineNow
                fgetl(fID);
                % Skip enough lines to get back to start of data
            end

            if nTabs == 3
                % is ZEN output file

                coordRead = textscan(fID, '%s\t%s\t%f\t%f');
                fclose(fID);

                obj.ROICoordinates = cell(1);

                    [~, IDstring, ~] = fileparts(obj.fileName);
                    % coordinates here are in "reslution units", which is ~10
                    % nm in most cases for ZEN output
                    thisCellsROIs = [coordRead{3}(strcmp(IDstring, cellstr(coordRead{2})))/scaleFactor, ...
                        obj.MaxSize - coordRead{4}(strcmp(IDstring, cellstr(coordRead{2})))/scaleFactor];

                    for p = 1:size(thisCellsROIs, 1)
                        obj.ROICoordinates{p} = zeros(5, 2);
                        % Assign ROI coordinates in proper format for inpolygon()
                        obj.ROICoordinates{p} = [thisCellsROIs(p,:) + [-obj.ROISize/2 -obj.ROISize/2];
                            thisCellsROIs(p,:) + [obj.ROISize/2 -obj.ROISize/2];
                            thisCellsROIs(p,:) + [obj.ROISize/2 obj.ROISize/2];
                            thisCellsROIs(p,:) + [-obj.ROISize/2 obj.ROISize/2];
                            thisCellsROIs(p,:) + [-obj.ROISize/2 -obj.ROISize/2]];
                    end

                loadOK = true;

    elseif nTabs > 3
        % Is polygonal format file
        % Each line may have a different number of coordinates, but
        % should always be paired
        % Everything is still in "resolution units", so be sure to
        % incorporate scaleFactor into the import

        roiCount = 1;
        fileEnd = false;
        cellList = cell(1,1);
        coordList = cell(1,1);
        while ~fileEnd
            thisLine = fgetl(fID);
            if ischar(thisLine)
                thisLine = strsplit(thisLine, sprintf('\t'));
                cellList{roiCount} = thisLine{2};
                coordList{roiCount} = reshape(str2double(thisLine(3:end)), 2, [])';
                roiCount = roiCount + 1;
            else
                fileEnd = true;
            end
        end

        fclose(fID);
        
        obj.ROICoordinates = cell(1);

            [~, IDstring, ~] = fileparts(obj.fileName);
           
            thisCellsROIs = coordList(strcmp(IDstring, cellList));

%             disp(length(thisCellsROIs));
            
            for p = 1:length(thisCellsROIs)

                % Assign ROI coordinates in proper format for inpolygon()
                obj.ROICoordinates{p} = thisCellsROIs{p}([1:end 1], :)/scaleFactor;
                obj.ROICoordinates{p}(:,2) = obj.MaxSize - obj.ROICoordinates{p}(:,2);

            end

        loadOK = true;                
                
            else
                error('Import file format not supported');
            end
        end
%-------------------------------------------------------------------------%   
        % Load ROI coordinates from coordinates.txt file (if existing)
        function v = Define_Square_ROIs_Auto(obj,~) % side of a square, um
            
                    anm = obj.Square_ROIs_Auto_anm;
                    qthresh = obj.Square_ROIs_Auto_qthresh;
                    Nrois = obj.Square_ROIs_Auto_maxNrois;
            
            downscale_factor = obj.pixelSizenm; % just widefield
            %downscale_factor = obj.pixelSizenm/10; % to see better use /10);                        
            XcT = obj.CellData(:,5);
            YcT = obj.CellData(:,6);
               YcT = max(YcT)-YcT;             
            imW = obj.MaxSize;
            imH = obj.MaxSize;
                SY = ceil(imW/downscale_factor);
                SX = ceil(imH/downscale_factor);   
                ash = ASH_2d(SX,SY,[YcT XcT]/downscale_factor,3);
                v = ash;                
            % to have pixel roughly the size of ROI    
            f = anm/obj.pixelSizenm;
            z = imresize(v,1/f);          
            z = z.*(z > quantile(z(:),qthresh));
            v = z;
            %
            obj.ROICoordinates = cell(0);
            for k=1:size(z,1)
                for m=1:size(z,2)
                    if z(k,m)>0
                       % disp([k-0.5,m-0.5]);
                      y = round(( (size(z,1)-k+1) -.5)*f*obj.pixelSizenm);
                      x = round((m-.5)*f*obj.pixelSizenm);
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
            % if too many ROIs, choose random Nrois among defined
            if (numel(obj.ROICoordinates) > Nrois)
                obj.ROICoordinates = obj.ROICoordinates(randi(numel(obj.ROICoordinates),1,Nrois));
            end
            %
        end
   %-------------------------------------------------------------------------%        
    end % methods
            
end
