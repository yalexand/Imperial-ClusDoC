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
    
    % ICA (intensity correlation) - the idea is from
        % J Neurosci. 2004 Apr 21;24(16):4070-81. 
        % A Syntaxin 1, Galpha(o), and N-type Calcium Channel Complex at a Presynaptic Nerve Terminal: Analysis by Quantitative Immunocolocalization
        % Qi Li 1, Anthony Lau, et al   
    mask = n1>0 & n2>0;
    corrimg = n1.*n2.*mask;
    s = corrimg(mask);
    ICA_score = mean(log(s(:)));
    %
    v = zeros(size(n1,1),size(n1,2),3,1,1);
    v(:,:,1,1,1) = n1;
    v(:,:,2,1,1) = n2;
    v(:,:,3,1,1) = log(1+corrimg);
        
    max_size = 49; % 7x7 pixels ~ 350x350 nm;
    n1_bigs_mask = bwareaopen(n1>0,max_size);
    n2_bigs_mask = bwareaopen(n2>0,max_size);
                          
    info = ['JUST jaccard : ' num2str(jaccard(n1>0,n2>0)), ... 
                 '; BIGS jaccard : ' num2str(jaccard(n1_bigs_mask,n2_bigs_mask)), ...
                 '; ICA score : ' num2str(ICA_score), ...
                 ];

end

