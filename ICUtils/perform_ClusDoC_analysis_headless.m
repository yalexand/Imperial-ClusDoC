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

ac.Define_Square_ROIs_Auto;

ac.Analyze_ROIs_DBSCAN(true); % verbose
ac.Analyze_ROIs_RipleyK; 
ac.Analyze_ROIs_DoC;

disp(['execution time ' num2str(toc/60) ' min']);

end


