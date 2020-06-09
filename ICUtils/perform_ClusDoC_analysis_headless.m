function perform_ClusDoC_analysis_headless(src_dir,filename,dst_dir,full_settings_file_name)

tic

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
ac.Load_Data(filename,src_dir);
ac.Define_Square_ROIs_Auto;
ac.Analyze_ROIs_DBSCAN(true); % verbose
ac.Analyze_ROIs_RipleyK; 

disp(['execution time ' num2str(toc/60) ' min']);

end


