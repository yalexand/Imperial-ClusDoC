%-------------------------------------------------------------%
function resolve_xlwrite_path_issue()
            jPath = javaclasspath('-dynamic');
            % first check it isn't already in the dynamic path
            WriteXLPath = false;
            for i = 1:length(jPath)
                if strfind(jPath{i},'jxl.jar');
                    WriteXLPath = true;
                    break;
                end
            end
                
            if ~WriteXLPath
                path = which('jxl.jar');
                if isempty(path)
                    path = fullfile(fileparts(mfilename('fullpath')), 'jxl.jar');
                end
                if ~isempty(path) && exist(path, 'file') == 2
                    javaaddpath(path);
                else 
                     assert('Cannot automatically locate an jxl JAR file');
                end
            end
            
             % first check it isn't already in the dynamic path
            WriteXLPath = false;
            for i = 1:length(jPath)
                if strfind(jPath{i},'MXL.jar');
                    WriteXLPath = true;
                    break;
                end
            end
                
            if ~WriteXLPath
                path = which('MXL.jar');
                if isempty(path)
                    path = fullfile(fileparts(mfilename('fullpath')), 'MXL.jar');
                end
                if ~isempty(path) && exist(path, 'file') == 2
                    javaaddpath(path);
                else 
                     assert('Cannot automatically locate an MXL JAR file');
                end
            end           
end