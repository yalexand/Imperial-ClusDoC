
classdef ClusDoC_MIiSR_analysis_controller < ClusDoC_analysis_controller 
    
    % 2021 Imperial College London.
    properties(SetObservable = true)            
        
        % these structures were copied from MIiSR interface, there are lot of
        % excessive stuff..
        MIiSR_Specs = struct(...
               'path','C:\', ...
                'nCh',2, ...
            'Ch1Name','Name for Channel 1', ...
            'Ch1Path','C:\Users\alexany\Imperial-ClusDoC\sequence_ThunderSTROM_May_22_2019_ch1.mat', ...
            'Ch2Name','Name for Channel 2', ...
            'Ch2Path','C:\Users\alexany\Imperial-ClusDoC\sequence_ThunderSTROM_May_22_2019_ch2.mat', ...
            'Ch3Name',[], ...
            'Ch3Path',[], ...
            'fMode',[1 1], ...
      'densityFilter',[false false], ...
         'densitySD',[NaN NaN], ...      
         'densitySD1',NaN, ...
         'densitySD2',NaN, ...
         'densitySD3',NaN, ...
        'densityDist',[200. 200.], ...
           'SAAcheck',1, ...
            'SAAdist','500', ...
            'SAArand','5', ...
              'SAAaF','1', ...
           'polyRand',0, ...
            'calcCDC',1, ...
             'CDC_SD','1.65', ...
          'SAAcutoff','20', ...
         'SAAuserCDC',0, ...
    'SAAuserCDCvalue','30', ...
           'RDFcheck',1, ...
        'spatialDist','500', ...
         'spatialCh1',1, ...
         'spatialCh2',1, ...
        'DBSCANcheck',1, ...
            'DBSCANk',[12 12], ...
            'DBSCANe',[12 12], ...
           'DBSCANch',1, ...
        'OPTICScheck',1, ...
            'OPTICSk',[12 12] , ...
           'OPTICSch',1, ...
         'OPTICShier',1, ...
          'hierRatio','0.75', ...
         'cropCoords',[149.2120 492.6886 185.3327 636.2657], ...
              'fName',{'Name for Channel 1_Name for Channel 2 - [149.212, 492.6886, 185.3327, 636.2657]'} ...        
              );
          
        MIiSR_Conditions = struct(...
              'fName','Name for Channel 1_Name for Channel 2 ', ...
             'ImScale',20, ...
            'saveTIFF',1, ...
       'saveUncropped',0, ...
          'minPhotons',0, ...
       'minPrecission',0, ...
        'savePosFiles',1, ...
        'saveSingleCh',0, ...
            'Ch1Label','Name for Channel 1', ...
            'Ch2Label','Name for Channel 2', ...
            'Ch3Label','', ...
                  'SD',1.6500, ...
              'Cutoff',20, ...
              'BinMax',500, ...
          'Iterations',5, ...
    'AnalyzedFraction',1, ...
            'polyRand',0, ...
             'gFilter',3, ...
      'ContourDensity',25, ...
               'Bound',500, ...
        'spatialGraph',1 ...
        );
        
    NA = 0.7;
    lambda = [ 647 555 0]; %channels 
    GL2photons = 1;    
    
    MIiSR_data = cell(2,1);
    
    end % properties
              
    methods

        function obj = ClusDoC_MIiSR_analysis_controller(varargin)
           obj = obj@ClusDoC_analysis_controller(varargin);
        end
                        
        %-------------------------- provide input data in MIiSR format
        function Load_Data(obj,fileName,pathName,channel)
            
                Load_Data@ClusDoC_analysis_controller(obj,fileName,pathName,channel);
                                        
                U_ind = -1; 
                if strcmp(char(obj.Original_from_file_header{channel}{1}(3)),'"x [pix]"') % WindSTORM
                    x_ind = 3;
                    y_ind = 4;
                    I_ind = 5;
                    spatial_multiplier = 1/obj.pixelSizenm; % corrections are in nanometers but WindSTORM format is in pixels
                elseif strcmp(char(obj.Original_from_file_header{channel}{1}(2)),'"x [nm]"') % ThunderSTORM
                    x_ind = 2;
                    y_ind = 3;
                    I_ind = 5;
                    U_ind = 8;                    
                    spatial_multiplier = 1;            
                elseif strcmp(char(obj.Original_from_file_header{channel}{1}(3)),'"x [nm]"') % another ThunderSTORM
                    x_ind = 3;
                    y_ind = 4;
                    I_ind = 6;
                    U_ind = 10;                    
                    spatial_multiplier = 1;
                else % X3
                end

                x = obj.Original_from_file_data{channel}{1,x_ind}*spatial_multiplier;
                y = obj.Original_from_file_data{channel}{1,y_ind}*spatial_multiplier;
                z = zeros(size(obj.Original_from_file_data{channel}{1,3}));
                I = obj.Original_from_file_data{channel}{1,I_ind}*obj.GL2photons;
                %
                if -1~=U_ind
                    uncertainty = obj.Original_from_file_data{channel}{1,U_ind};
                else % calculate via PSF width and Intensity
                    PSFwidth = SMLM_sigma_calculator(obj.NA,obj.lambda(channel),obj.pixelSizenm)*obj.pixelSizenm;
                    uncertainty = PSFwidth./sqrt(I);
                end
                %
                % one of ThnderSTORM fomats may be not corrected for that
                % under loading
                y = max(y) - y;
                %
                obj.MIiSR_data{channel} = [ x y z I uncertainty];                                                
        end
        
        %-------------------------- apply registration to NIiSR data
        function Apply_channel2_registration_corrections(obj,dx2dy2,~)            
            Apply_channel2_registration_corrections@ClusDoC_analysis_controller(obj,dx2dy2);            
            X = obj.MIiSR_data{2}(:,1) - dx2dy2(1);
            Y = obj.MIiSR_data{2}(:,2) - dx2dy2(2);
            nmppix = obj.pixelSizenm;
            mask  = X>=0 & X<=obj.SizeX*nmppix & Y>=0 & Y<=obj.SizeY*nmppix;
            obj.MIiSR_data{2} = obj.MIiSR_data{2}(mask,:); 
            obj.MIiSR_data{2}(:,1) = X(mask);
            obj.MIiSR_data{2}(:,2) = Y(mask);                        
        end
        
        %-------------------------- 
        function ROI_data_MIiSR = get_ROI_data_MIiSR(obj,roi_index,channel)
            ROI_data_MIiSR = [];
            try
                roi = obj.ROICoordinates{roi_index};                        
                x_roi=roi(:,1);
                y_roi=roi(:,2);
                x = obj.MIiSR_data{channel}(:,1);
                y = obj.MIiSR_data{channel}(:,2);
                whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi); 
                ROI_data_MIiSR = obj.MIiSR_data{channel}(whichPointsInROI,:);
                %
                % MIiSR convention
                ROI_data_MIiSR(:,1) = ROI_data_MIiSR(:,1)-min(ROI_data_MIiSR(:,1));
                ROI_data_MIiSR(:,2) = ROI_data_MIiSR(:,2)-min(ROI_data_MIiSR(:,2));                
            catch
                disp('get_ROI_data_MIiSR: error!');
            end                        
        end
             
        %-------------------------- 
        function SAAstruct = SAA(obj,~)
               % 
               SAA_out_dirname = cell(2,1);
               for chan = 1 : 1 %obj.Nchannels                   
                   fname = strrep(obj.fileName{chan},'.csv','');
                   fname = strrep(fname,'.txt','');
                   SAA_out_dirname{chan} = [obj.Outputfolder filesep fname '_channel_' num2str(chan) '_MIiSR_Results' filesep  'SAA'];
                   if ~exist(SAA_out_dirname{chan},'dir')
                       mkdir( fullfile(obj.Outputfolder,[fname '_channel_' num2str(chan) '_MIiSR_Results'],'SAA'));
                   end
               end
               %
               SAAstruct = cell(length(obj.ROICoordinates),1);
               for roiIter = 1:length(obj.ROICoordinates) % ROI number                        
                    disp([roiIter length(obj.ROICoordinates)]);
                    try
                        SAAstruct{roiIter} = obj.SAA2col(roiIter,SAA_out_dirname{chan});
                    catch
                        disp('SAA2col - error');
                        SAAstruct{roiIter} = [];
                    end
               end
               save([SAA_out_dirname{chan} filesep 'SAAstruct'],'SAAstruct');
        end
        
        %-------------------------- 
        function SpatialData = SpatialStatistics(obj,~)
               % 
               SpD_out_dirname = cell(2,1);
               for chan = 1 : obj.Nchannels
                   fname = strrep(obj.fileName{chan},'.csv','');
                   fname = strrep(fname,'.txt','');
                   SpD_out_dirname{chan} = [obj.Outputfolder filesep fname '_channel_' num2str(chan) '_MIiSR_Results' filesep  'SpatialStats'];
                   if ~exist(SpD_out_dirname{chan},'dir')
                       mkdir( fullfile(obj.Outputfolder,[fname '_channel_' num2str(chan) '_MIiSR_Results'],'SpatialStats'));
                   end
               end
               %
               for chan = 1 : obj.Nchannels
                   SpatialData = cell(length(obj.ROICoordinates),1);
                   for roiIter = 1:length(obj.ROICoordinates) % ROI number                        
                        disp([roiIter length(obj.ROICoordinates)]);
                        try
                            if 1==chan
                                Ch1 = 1;
                                Ch2 = 2;
                            else
                                Ch1 = 2;
                                Ch2 = 1;                                
                            end
                            [SpatialData{roiIter},~] = obj.spatialStats(roiIter,SpD_out_dirname{chan},Ch1,Ch2);
                        catch
                            disp('spatialStats - error');
                            SpatialData{roiIter} = [];
                        end
                   end
                   save([SpD_out_dirname{chan} filesep 'SpatialData'],'SpatialData');
               end
        end        
        
        %-------------------------- 
        function MIiSR_DBSCAN(obj,~)
               % 
               DBSCAN_out_dirname = cell(2,1);
               for chan = 1 : obj.Nchannels
                   fname = strrep(obj.fileName{chan},'.csv','');
                   fname = strrep(fname,'.txt','');
                   DBSCAN_out_dirname{chan} = [obj.Outputfolder filesep fname '_channel_' num2str(chan) '_MIiSR_Results' filesep  'DBSCANStats'];
                   if ~exist(DBSCAN_out_dirname{chan},'dir')
                       mkdir( fullfile(obj.Outputfolder,[fname '_channel_' num2str(chan) '_MIiSR_Results'],'DBSCANStats'));
                   end
               end
               %
               for chan = 1 : obj.Nchannels
                   DBSCANData = cell(length(obj.ROICoordinates),1);
                   for roiIter = 1:length(obj.ROICoordinates) % ROI number                        
                        disp([roiIter length(obj.ROICoordinates)]);
%                        try                     
                            [DBSCANout.DBSCANmat,DBSCANout.DBSCANtable,~] = obj.DBSCANStats(roiIter,chan,DBSCAN_out_dirname{chan});
                            DBSCANData{roiIter} = DBSCANout;
%                         catch
%                             disp('DBSCANStats - error');
%                             DBSCANData{roiIter} = [];
%                         end
                   end               
                   save([DBSCAN_out_dirname{chan} filesep 'DBSCANData'],'DBSCANData');
               end            
        end

        %-------------------------- 
        function MIiSR_OPTICS(obj,~)
               % 
               OPTICS_out_dirname = cell(2,1);
               for chan = 1 : obj.Nchannels
                   fname = strrep(obj.fileName{chan},'.csv','');
                   fname = strrep(fname,'.txt','');
                   OPTICS_out_dirname{chan} = [obj.Outputfolder filesep fname '_channel_' num2str(chan) '_MIiSR_Results' filesep  'OPTICSStats'];
                   if ~exist(OPTICS_out_dirname{chan},'dir')
                       mkdir( fullfile(obj.Outputfolder,[fname '_channel_' num2str(chan) '_MIiSR_Results'],'OPTICSStats'));
                   end
               end
               %
               for chan = 1 : obj.Nchannels
                   OPTICSData = cell(length(obj.ROICoordinates),1);
                   for roiIter = 1:length(obj.ROICoordinates) % ROI number                        
                        disp([roiIter length(obj.ROICoordinates)]);
                        try                                                 
                            [OPTICStable,Epsilon] = obj.OPTICSStats(roiIter,chan,OPTICS_out_dirname{chan});
                            OPTICSout{1} = OPTICStable;
                            OPTICSout{2} = Epsilon;
                            OPTICSData{roiIter} = OPTICSout;
                        catch
                            disp('OPTICSStats - error');
                            OPTICSData{roiIter} = [];
                        end
                   end               
                   save([OPTICS_out_dirname{chan} filesep 'OPTICSData'],'OPTICSData');
               end            
        end
                                      
    end % methods
end



















