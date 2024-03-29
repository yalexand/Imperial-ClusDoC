% Import1File will pull the data from a Zeiss 1.txt file given in string to
% filepath fname.  Output is structure with fields .Header, .Data, and
% .Footer.
% 
% 05/14/2017 added support for Nikon SMLM file import.  Data is extracted
% and placed into format to match 13-column Zeiss data file.

function [Data,Head_text,Body_text] = ImportFile(fname,handles)

filepath = fname;

fid = fopen(filepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format-specific import:
% First entry in data is a '1' in the first column.  This should be the
% first numeral following a whole string of characters and white space.

% Issue with new file formats: a 13th column of data (Z Slice) is included
% in some files.  Check if it is present by counting tabs in the first
% line.
[~, ~, ext] = fileparts(filepath);
if strcmp(ext, '.txt')
    TABCHAR = sprintf('\t');
elseif strcmp(ext, '.csv')
    TABCHAR = sprintf(',');
else
    error('File name extension not supported. File must be .txt or .csv');
end

idx = find(fgetl(fid) == TABCHAR);
fseek(fid, 0, -1);

% Pull headers
format_spec = '%s';
N_cols = numel(idx)+1; % Number of tabs in first line, plus 1

[Head_text, Head_post] = textscan(fid, format_spec, N_cols, 'delimiter', TABCHAR);

if (length(Head_text{1}) == 24) ||  (length(Head_text{1}) == 26) % Nikon file format
    format_spec = ['%s', repmat('%n', 1, N_cols - 1)]; 
else % Zeiss file format
    Body_format = '%n';
    format_spec = repmat(Body_format, 1, numel(idx)+1);
end
[Body_text, Body_post] = textscan(fid, format_spec, 'delimiter', TABCHAR);

% Odd end-of-file bit in Zeiss files I can't seem to replicate in
% MATLAB-generated files.  Going to just search for the first string in the
% footer and start from there. 
floc = [];
while 1
    ftloc = ftell(fid); % Get byte position of line you're on
    tline = fgetl(fid); % Read string from line
    if ~ischar(tline)
       break
    end
    
    if strfind(tline, 'VoxelSizeX');
        
        floc = strfind(tline, 'VoxelSizeX'); % Pull Start position of first string in footer
        
        break
    else
        floc = [];
    end
    
end

if ~isempty(floc)

    format_spec = '%s%f%s';
    % Try-catch block to solve file import error on some (old?) files
    try
        Footer_text = textscan(fid, format_spec, 'headerLines', 1, 'delimiter', ':');
        fseek(fid, (ftloc-floc), -1); % Footer starts at start of line where string match was found.

        VoxelSizeX = Footer_text{2}(1);
        VoxelSizeY = Footer_text{2}(2);
        ResolutionX = Footer_text{2}(3);
        ResolutionY = Footer_text{2}(4);
        SizeX = Footer_text{2}(5);
        SizeY = Footer_text{2}(6);
    
    catch
        % On some machines the 'headerLines' value has to be 0 to work.  Not sure
        % why this is the case.
        fseek(fid, (ftloc), -1);
        Footer_text = textscan(fid, format_spec, 'headerLines', 0, 'delimiter', ':');
        Footer_text{1}(end) = [];
        
        VoxelSizeX = Footer_text{2}(1);
        VoxelSizeY = Footer_text{2}(2);
        ResolutionX = Footer_text{2}(3);
        ResolutionY = Footer_text{2}(4);
        SizeX = Footer_text{2}(5);
        SizeY = Footer_text{2}(6);
    end
else
    Footer_text = '';
end

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

    switch length(Body_text)

        case 3 % simulated data format
                    %dhead = '"x","y","index"

%                 Nmin = inf; % ?????
%                 for k=1:length(Body_text)
%                     if length(Body_text{k})<Nmin
%                         Nmin=length(Body_text{k});
%                     end
%                 end                
%                 for k=1:length(Body_text)
%                     Body_text{k} = Body_text{k}(1:Nmin);
%                 end              
            
                Data = [(1:length(Body_text{1}))', ... % index
                        ones(length(Body_text{1}), 1), ... % frame                    
                        ones(length(Body_text{1}), 1), ...
                        zeros(length(Body_text{1}), 1), ...
                        Body_text{1}, ...
                        Body_text{2}, ... %pixelSizenm*(max(Body_text{4})-Body_text{4}), ...                        
                        zeros(length(Body_text{1}), 1), ... % precision (uncertainty_xy)
                        zeros(length(Body_text{1}), 1), ... % intensity
                        zeros(length(Body_text{1}), 1), ... % bckgstd
                        ones(length(Body_text{1}), 1), ... % chi square
                        zeros(length(Body_text{1}), 1), ... % psf width [nm]
                        ones(length(Body_text{1}), 1), ...
                        ones(length(Body_text{1}), 1)];
                                                            
            % simulated data format

        % YA - another ThunderSTORM format
        case {8,9} % ThunderSTORM format.  The user should have chosen the .csv file.  There's an associated XML-ish
                % .txt file that has metadata and at least the pixel size in it. 
                % This file must be the same name as the input but with
                % -protocol.txt appended.
                %
                shift = 0;
                if 9==length(Body_text)
                    shift = 1;
                end
                
                Nmin = inf; % ?????
                for k=1:length(Body_text)
                    if length(Body_text{k})<Nmin
                        Nmin=length(Body_text{k});
                    end
                end                
                for k=1:length(Body_text)
                    Body_text{k} = Body_text{k}(1:Nmin);
                end                
                %
                Data = [(1:Nmin)', ...
                        Body_text{1+shift}, ...                    
                        ones(Nmin, 1), ...
                        zeros(Nmin, 1), ...
                        Body_text{2+shift}, ... %                          
                        handles.pixelSizenm*handles.SizeY - Body_text{3+shift}, ... %max(Body_text{3})-Body_text{3}, ... % Body_text{3}, ...
                        Body_text{8+shift}, ... % precision (uncertainty_xy)
                        Body_text{5+shift}, ... % intensity
                        Body_text{7+shift}, ... % bckgstd
                        ones(Nmin, 1), ... % chi square
                        Body_text{4+shift}, ... % psf width (sigmas)
                        ones(Nmin, 1), ...
                        ones(Nmin, 1)];                
                                                                
%                 Data = [(1:length(Body_text{1}))', ...
%                         Body_text{1+shift}, ...                    
%                         ones(length(Body_text{1}), 1), ...
%                         zeros(length(Body_text{1}), 1), ...
%                         Body_text{2+shift}, ... %                          
%                         handles.pixelSizenm*handles.SizeY - Body_text{3+shift}, ... %max(Body_text{3})-Body_text{3}, ... % Body_text{3}, ...
%                         Body_text{8+shift}, ... % precision (uncertainty_xy)
%                         Body_text{5+shift}, ... % intensity
%                         Body_text{7+shift}, ... % bckgstd
%                         ones(length(Body_text{1}), 1), ... % chi square
%                         Body_text{4+shift}, ... % psf width (sigmas)
%                         ones(length(Body_text{1}), 1), ...
%                         ones(length(Body_text{1}), 1)];
                
                [folderPath, metaDataFile, ~] = fileparts(filepath);
                mdID = fopen(fullfile(folderPath, strcat(metaDataFile, '-protocol.txt')));
                
                if mdID ~= -1
                    % proceed
                    % Pull out pixel size
                    md = textscan(mdID, '%s', 'delimiter', '\n');
                    md = md{1};
                    pixMD = ~cellfun(@isempty, strfind(md, 'pixel'));
                    if sum(pixMD) ~= 1
                        error('Format of protocol file incorrect');
                    end

                    splitMD = textscan(md{pixMD}, '%s %f,');
                    pixelSizenm = splitMD{2}; 
                    handles.pixelSizenm = pixelSizenm;
                else
                    disp('Associated ThunderSTORM *-protocol.txt file not found. Using default for pixelSizenm');
                end
                % Example data was only single-channel data.  Assuming this
                % is true for all ThunderSTORM data from here and filling
                % in dummy '1' values for all channels.                 
                            
        case 5 % WindSTORM output as arranged in Imperial Biophotonics            
                    %dhead = '"id","frame","x [pix]","y [pix]","intensity [photon]"';        

            pixelSizenm = handles.pixelSizenm;
            psf_width = handles.WindSTORM_Sigmapix*handles.pixelSizenm;

                Nmin = inf; % ?????
                for k=1:length(Body_text)
                    if length(Body_text{k})<Nmin
                        Nmin=length(Body_text{k});
                    end
                end                
                for k=1:length(Body_text)
                    Body_text{k} = Body_text{k}(1:Nmin);
                end              
            
                Data = [(1:length(Body_text{1}))', ... % index
                        Body_text{2}, ... % frame                    
                        ones(length(Body_text{1}), 1), ...
                        zeros(length(Body_text{1}), 1), ...
                        pixelSizenm*Body_text{3}, ...
                        pixelSizenm*(handles.SizeY-Body_text{4}), ... %pixelSizenm*(max(Body_text{4})-Body_text{4}), ...                        
                        zeros(length(Body_text{1}), 1), ... % precision (uncertainty_xy)
                        Body_text{5}, ... % intensity
                        zeros(length(Body_text{1}), 1), ... % bckgstd
                        ones(length(Body_text{1}), 1), ... % chi square
                        psf_width*ones(length(Body_text{1}), 1), ... % psf width [nm]
                        ones(length(Body_text{1}), 1), ...
                        ones(length(Body_text{1}), 1)];
                                                            
            % WindSTORM format - ends
            
        % YA - ends            
        
        case {10, 11} % ThunderSTORM format.  The user should have chosen the .csv file.  There's an associated XML-ish
                % .txt file that has metadata and at least the pixel size in it. 
                % This file must be the same name as the input but with
                % -protocol.txt appended.
                
                Data = [Body_text{1}, Body_text{2}, ones(length(Body_text{1}), 1), ...
                        zeros(length(Body_text{1}), 1), Body_text{3}, Body_text{4}, ...
                        Body_text{10}, Body_text{6}, Body_text{8}, Body_text{9}, Body_text{5}, ...
                        ones(length(Body_text{1}), 1), ones(length(Body_text{1}), 1)];
                
                [folderPath, metaDataFile, ~] = fileparts(filepath);
                mdID = fopen(fullfile(folderPath, strcat(metaDataFile, '-protocol.txt')));
                
                if mdID ~= -1
                    try
                    % proceed
                    % Pull out pixel size
                    md = textscan(mdID, '%s', 'delimiter', '\n');
                    md = md{1};
                    pixMD = ~cellfun(@isempty, strfind(md, 'pixel'));
                    if sum(pixMD) ~= 1
                        error('Format of protocol file incorrect');
                    end

                    splitMD = textscan(md{pixMD}, '%s %f,');
                    pixelSizenm = splitMD{2};
                    handles.pixelSizenm = pixelSizenm;
                    catch err
                        disp('cannot open protocol file format');
                        disp(err.message);
                    end
                else
                    disp('Associated ThunderSTORM *-protocol.txt file not found. Using default for pixelSizenm');
                end
                                            
        case 12

            Data = [Body_text{1} Body_text{2} Body_text{3} Body_text{4} Body_text{5} Body_text{6}...
                Body_text{7} Body_text{8} Body_text{9} Body_text{10} Body_text{11} Body_text{12}];

        case 13

            Data = [Body_text{1} Body_text{2} Body_text{3} Body_text{4} Body_text{5} Body_text{6}...
                Body_text{7} Body_text{8} Body_text{9} Body_text{10} Body_text{11} Body_text{12} Body_text{13}];
            Data(:,12) = ones(size(Data(:,12))); % YA: otherwise, problems

        case 14
            
            Data = [Body_text{1} Body_text{2} Body_text{3} Body_text{4} Body_text{5} Body_text{6}...
                Body_text{7} Body_text{8} Body_text{9} Body_text{10} Body_text{11} Body_text{12} Body_text{13} Body_text{14}]; 
            
        case 15 % Files generated by Generate3DTestData give 15 columns for some reason

            Data = [Body_text{1} Body_text{2} Body_text{3} Body_text{4} Body_text{5} Body_text{6}...
                Body_text{7} Body_text{8} Body_text{9} Body_text{10} Body_text{11} Body_text{12} Body_text{13} Body_text{14}]; 
            
        case {24, 26} % Nikon format

            % Need to get channel ID out of a string ID in the first
            % column
            
            % If 26 columns, Zw and Zwc are ignored

            channelID = ones(length(Body_text{1}), 1);
            possibleChannels = unique(Body_text{1});
            possibleChannels(cellfun(@isempty, possibleChannels)) = [];
            if length(unique(Body_text{1})) > 1
                for k = 1:length(possibleChannels)
                    channelID(strcmp(Body_text{1}, possibleChannels{k}), 1) = k;
                end
            end
            
            Data = [(1:(length(Body_text{1})))', Body_text{13}, Body_text{14}, zeros(length(Body_text{1}), 1), ...
                    Body_text{23}, Body_text{24}, Body_text{19}, Body_text{20}, Body_text{11}, zeros(length(Body_text{1}), 1), ...
                    Body_text{8}, channelID, Body_text{18}];                
    end

end