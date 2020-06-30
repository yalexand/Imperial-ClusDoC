function test_perform_ClusDoC_analysis_headless()

% SINGLE CHANNEL
% 	src_dir = ['.' filesep 'TestData'];
% 	dst_dir = src_dir; % :)    
%     filename = 'sequence_ThunderSTROM_May_22_2019_ch1.csv';
%     full_settings_file_name = [src_dir filesep 'sequence_ThunderSTROM_May_22_2019_settings.xml'];    
% 	perform_ClusDoC_analysis_headless(dst_dir,full_settings_file_name,src_dir,filename);
 
% TWO CHANNELS
	src_dir_ch1 = ['.' filesep 'TestData'];
    filename_ch1 = 'sequence_ThunderSTROM_May_22_2019_ch1.csv';
	src_dir_ch2 = ['.' filesep 'TestData'];
    filename_ch2 = 'sequence_ThunderSTROM_May_22_2019_ch2.csv';
        %
        dst_dir = src_dir_ch1; % :)
        full_settings_file_name = [src_dir_ch1 filesep 'sequence_ThunderSTROM_May_22_2019_settings.xml'];
        %
    perform_ClusDoC_analysis_headless(dst_dir,full_settings_file_name, ... 
        src_dir_ch1,filename_ch1, ...
        src_dir_ch2,filename_ch2);
    
%  % TWO CHANNELS - WindSTORM
% 	src_dir_ch1 = ['.' filesep 'TestData'];
%     filename_ch1 = 'Tubulin_sequence_WindSTORM_ch1.csv';
% 	src_dir_ch2 = ['.' filesep 'TestData'];
%     filename_ch2 = 'Tubulin_sequence_WindSTORM_ch2.csv';
%         %
%         dst_dir = src_dir_ch1; % :)
%         full_settings_file_name = [src_dir_ch1 filesep 'Tubulin_sequence_WindSTORM_settings.xml'];
%         %
%     perform_ClusDoC_analysis_headless(dst_dir,full_settings_file_name, ... 
%         src_dir_ch1,filename_ch1, ...
%         src_dir_ch2,filename_ch2);   
    

end