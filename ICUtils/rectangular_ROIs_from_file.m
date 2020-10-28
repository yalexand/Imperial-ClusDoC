function rois = square_ROIs_from_file(fullfilename)
% either "csv" (FIJI), or "xls" (ICY)
% output format: X0,Y0,Width,Height

    rois = [];

    [~,~,fext] = fileparts(fullfilename);
   
    if contains(fext,'csv') % FIJI

        d = importdata(fullfilename,',');
        cells = d.textdata;
        for ind=1:size(cells,2)
            if strcmp('X',cells(1,ind)), break, end
        end
        rois = zeros(size(cells,1)-1,4);
            for r=1:4 % row
                for c=1:size(cells,1)-1 % column
                  rois(c,r) = str2double(cells(c+1,ind+r-1));
                end
            end        
    elseif contains(fext,'xls') % ICY
            T = readtable(fullfilename);
            colnames = T.Properties.VariableNames;
            for ind=1:length(colnames)
                if strcmp('PositionX',colnames{ind}), break, end
            end
            
            rois = zeros(size(T,1),4);
            
                for c=1:size(T,1)% column
                  rois(c,1) = round(str2double(T{c,ind}));
                  rois(c,2) = round(str2double(T{c,ind+1}));
                  rois(c,3) = round(str2double(T{c,ind+5}));
                  rois(c,4) = round(str2double(T{c,ind+6}));
                end
    else
        disp('ROIs empty - file format not supported, load either "csv" (FIJI), or "xls" (ICY)');
    end                                            
end