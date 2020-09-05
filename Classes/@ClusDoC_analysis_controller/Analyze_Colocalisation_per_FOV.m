function [v,info,out_data] = Analyze_Colocalisation_per_FOV(obj,~)

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
    %ICA_score = mean(log(s));
    ICA_score = mean(log(s))/std(log(s)); % this works as well
    %ICA_score = mean(s)/std(s);          % this doesn't work (?) 
    corrimg = log(1+corrimg);

    % this gives nice ICA image
%     s1 = n1(mask);
%     s2 = n2(mask);
%     corrimg = (n1-mean(s1)).*(n2-mean(s2))/std(s1)/std(s2);
%     corrimg(mask==0)=0;
%     % use multiplicative way on corrimg>0 to get co-localization quantifier  
%     s = n1.*n2.*(corrimg>0);
%     s = s(s~=0);
%     ICA_score = mean(log(s))/std(log(s));  
    %
    v = zeros(size(n1,1),size(n1,2),3,1,1);
    v(:,:,1,1,1) = n1;
    v(:,:,2,1,1) = n2;
    v(:,:,3,1,1) = corrimg;
        
    max_size = 49; % 7x7 pixels ~ 350x350 nm;
    n1_bigs_mask = bwareaopen(n1>0,max_size);
    n2_bigs_mask = bwareaopen(n2>0,max_size);
                          
    info = ['JUST jaccard : ' num2str(jaccard(n1>0,n2>0)), ... 
                 '; BIGS jaccard : ' num2str(jaccard(n1_bigs_mask,n2_bigs_mask)), ...
                 '; ICA score : ' num2str(ICA_score), ...
                 ];
    out_data = [jaccard(n1>0,n2>0) jaccard(n1_bigs_mask,n2_bigs_mask) ICA_score];
end

