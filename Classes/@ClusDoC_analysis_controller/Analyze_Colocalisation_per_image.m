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
    
    mask = n1>0 | n2 >0;
    s1 = n1(mask);
    s2 = n2(mask);
    corr_or = corr(s1(:),s2(:));
   
    % ICA (intensity correlation)
        % J Neurosci. 2004 Apr 21;24(16):4070-81. 
        % A Syntaxin 1, Galpha(o), and N-type Calcium Channel Complex at a Presynaptic Nerve Terminal: Analysis by Quantitative Immunocolocalization
        % Qi Li 1, Anthony Lau, et al      
    z1 = map(n1,0,1);
    z2 = map(n2,0,1);
        s1 = z1(mask);
        s2 = z2(mask);
            t1 = quantile(s1(:),.99);
            t2 = quantile(s2(:),.99);
                z1(z1>t1)=t1;
                z2(z2>t2)=t2;
                    s1 = z1(mask);
                    s2 = z2(mask);
                        m1 = mean(s1(:));
                        m2 = mean(s2(:));
                        corrimg = (z1-m1).*(z2-m2);
                        t = mean(corrimg(~mask));
                        corrimg(corrimg<t)=0;
                        s = corrimg(0~=corrimg);
                        ICA_score = mean(s(:));
    %
    v = zeros(size(n1,1),size(n1,2),3,1,1);
    v(:,:,1,1,1) = n1;
    v(:,:,2,1,1) = n2;
    v(:,:,3,1,1) = corrimg;
        
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
                 '; ICA score : ' num2str(ICA_score), ...
                 ];

end

