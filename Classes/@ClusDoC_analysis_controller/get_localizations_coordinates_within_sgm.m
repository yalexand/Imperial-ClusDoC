function data = get_localizations_coordinates_within_sgm(obj,chan,mag,~)
%Localizations_Coordinates_within_sgm 
%   mag - magnification between 1nm/pixel, and segmentation mask
%   obj.sgm is presumed binary    
    data = [obj.CellData{chan}(:,5), obj.CellData{chan}(:,6)];
    data_red = round(data/mag);
    N = size(data,1);
    mask = ones(N,1);
    for k=1:N
        x = data_red(k,1);
        y = data_red(k,2);
        %
        if x<1 || x>obj.SizeX || y<1 || y>obj.SizeY || 0==obj.sgm(y,x)
            mask(k)=0;
        end
    end
    data(mask==0,:)=[];
end

