function [v,info] = Analyze_Colocalisation_per_image(obj,~)

    pixel_size_in_microns = 0.05; % 50 nanometers
    
    try            
        % localisations numbers per pixel
        n1 = obj.get_localisation_numbers(1,pixel_size_in_microns);
        n2 = obj.get_localisation_numbers(2,pixel_size_in_microns);
    catch err
        disp(err.message);
        return;
    end
    %

    %
    v = zeros(size(n1,1),size(n1,2),2,1,1);
    v(:,:,1,1,1) = n1;
    v(:,:,2,1,1) = n2;
    
    mask = n1>0 | n2 >0;
    s1 = n1(mask);
    s2 = n2(mask);
    corr_or = corr(s1(:),s2(:));

    max_size = 49; % 7x7 pixels ~ 350x350 nm;
    n1_bigs_mask = bwareaopen(n1>0,max_size);
    n2_bigs_mask = bwareaopen(n2>0,max_size);
    
    %icy_imshow(n1_bigs_mask + n2_bigs_mask);
    %icy_imshow(n1_bigs_mask & n2_bigs_mask);
    
%     bigs_or_mask = n1_bigs_mask(:) | n2_bigs_mask(:);
%                     
%     corr_or_bigs= corr(s1(bigs_or_mask),s2(bigs_or_mask));       
%     
%     info = ['OR : ' num2str(corr_or), ...
%                  '; AND : ' num2str(corr_and), ...
%                  '; BIGS : ' num2str(corr_or_bigs), ... 
%                  '; BIGS jaccard : ' num2str(jaccard(n1_bigs_mask,n2_bigs_mask))];

    bigs_or_mask = n1_bigs_mask(:) | n2_bigs_mask(:);                        
    info = ['JUST jaccard : ' num2str(jaccard(n1>0,n2>0)), ... 
                 '; BIGS jaccard : ' num2str(jaccard(n1_bigs_mask,n2_bigs_mask)), ...
                 '; corr or : ' num2str(corr_or), ...
                 ];

end

