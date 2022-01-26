%-------------------------------------------------------------------------%
        function Load_simulated_CAML_tsv_data(obj,fileName,pathName,pixelSizenm,chan,~)
            
            obj.fileName{chan} = fileName;
            obj.pathName{chan} = pathName;
                
            if contains(obj.fileName{chan},'.tsv')
                
   % Assign variables pulled out of each step.
    % Output table needs to be in order:
    % 1.  Index	
    % 2.  First Frame	
    % 3.  Number Frames	
    % 4.  Frames Missing	
    % 5.  Position X [nm]	
    % 6.  Position Y [nm]	
    % 7.  Precision [nm]	
    % 8.  Number Photons	
    % 9.  Background variance	
    % 10. Chi square	
    % 11. PSF width [nm]	
    % 12. Channel	
    % 13. Z Slice   
                
            obj.pixelSizenm = pixelSizenm;
            psf_width = 150; % :)
            
            [data,original_header,original_data] = tsvread(fullfile(obj.pathName{chan},obj.fileName{chan}));

                N = size(data,1);    
                importData = [(1:N)', ... % index
                        ones(N,1), ... % frame                    
                        ones(N,1), ... % 
                        zeros(N,1), ...
                        data(:,1), ...
                        data(:,2), ...                        
                        zeros(N,1), ... % precision (uncertainty_xy)
                        zeros(N,1), ... % intensity
                        zeros(N,1), ... % bckgstd
                        ones(N,1), ... % chi square
                        psf_width*ones(N,1), ... % psf width [nm]
                        ones(N,1), ...
                        ones(N,1)];
                                                                                                                                                               
                obj.CellData{chan} = [importData zeros(size(importData, 1), 8)];

                obj.NDataColumns = size(importData, 2);
                
                obj.CellData{chan}(:,obj.NDataColumns + 2) = 1; % All data is in mask until set otherwise
                
% ?                
%                 obj.CellData{chan}(any(obj.CellData{chan}(:, 5) > obj.SizeX), : )= [];
%                 obj.CellData{chan}(any(obj.CellData{chan}(:, 6) > obj.SizeY), : )= [];                
%                 obj.CellData{chan}(any(obj.CellData{chan}(:, 5:6) < 0), : )= [];
                
                obj.Nchannels = max(chan,obj.Nchannels); %min([numel(unique(obj.CellData{chan}(:,12))), 2]); % cap import to 2 channels ever
                
                obj.Original_from_file_data{chan} = original_data;
                obj.Original_from_file_header{chan} = original_header;
                
            else                
                fprintf(1, 'File not in accepted coordinate table format.\nSkipping %s\n', fullfile(obj.pathName{chan},obj.fileName{chan}));
            end            
