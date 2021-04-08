function perform_ClusDoC_analysis_headless(dst_dir,full_settings_file_name,varargin)

tic

src_dir1 = varargin{1};
filename1 = varargin{2};

try
    src_dir2 = varargin{3};
    filename2 = varargin{4};
catch
end
    
ac = ClusDoC_analysis_controller();

if isfile(full_settings_file_name)
    try
    ac.load_settings(full_settings_file_name);
    catch ex
        disp(ex.message());
        disp('wrong or corrupt settings file, cannot continue');
        return;
    end
else
    disp('wrong settings file path, cannot continue');
    return;
end

ac.Outputfolder = dst_dir;   

ac.Load_Data(filename1,src_dir1,1);
if exist('filename2','var') && exist('src_dir2','var')
    ac.Load_Data(filename2,src_dir2,2);
    % ac.Align_channels_nmppix = 10;
    % ac.Align_channels_method = 'xcorr2';    
    % [dx2,dy2] = ac.Get_channel2_registration_corrections
    % ac.Apply_channel2_registration_corrections([dx2 dy2]);
end

% for Tubulin sequence_WindSTORM.csv - 256x256, 107nm, Sigma 1.2
%ac.Square_ROIs_Auto_anm = 3000;
%ac.Square_ROIs_Auto_qthresh = [.8 .8];

%ac.Square_ROIs_Auto_qthresh = [.7 .7];
%ac.Square_ROIs_Auto_method = 'composite';
%ac.Square_ROIs_Auto_verbose = false;
%ac.Define_Square_ROIs_Auto;

ac.Square_ROIs_Auto_anm = 2500;
ac.Square_ROIs_Auto_method = 'composite';
ac.Square_ROIs_Auto_verbose = true;
ac.Load_Segmentation([src_dir1 filesep 'sequence_ThunderSTROM_May_22_2019_sgm.tif']);
h = ac.Define_Square_ROIs_From_Segmentation; 

%ac.Analyze_ROIs_DBSCAN(true); % verbose or not

%ac.Analyze_ROIs_RipleyK; 

%ac.Analyze_ROIs_DoC;

% dx2 = 170;
% dy2 = -90;
% ac.Save_original_channel2_data_with_XY_registration_corrections(ac.Outputfolder,[dx2 dy2]);

disp(['execution time ' num2str(toc/60) ' min']);

end


