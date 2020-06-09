function test_perform_ClusDoC_analysis_headless()

	src_dir = ['.' filesep 'TestData'];
    
    filename = 'sequence_ThunderSTROM_May_22_2019.csv';

	dst_dir = src_dir; % :)

	full_settings_file_name = [src_dir filesep 'sequence_ThunderSTROM_May_22_2019_settings.xml'];

	perform_ClusDoC_analysis_headless(src_dir,filename,dst_dir,full_settings_file_name);

end