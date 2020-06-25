function tform = Correct_Simple_Image_Translation_xcorr2(moving,static,s1,s2)
%Correct_Simple_Image_Translation_xcorr2 
%   function operates on a combination of original and its 3-size gaussian derivatives      

        [gx,gy] = gsderiv(moving,s1,1);
            moving_comp = sqrt(gx.*gx + gy.*gy) + moving/3;
            
        [gx,gy] = gsderiv(static,s1,1);
            static_comp = sqrt(gx.*gx + gy.*gy) + static/3;

        z = xcorr2_fft(moving_comp,static_comp); 
        %
        g1 = gsderiv(z,s1,0);
        g2 = gsderiv(z,s2,0);
        z = (g1-g2);
        %
        z(z<0)=0;
        %
        [wc,hc] = size(z);
        wc=fix(wc/2);
        hc=fix(hc/2);
        rx = wc-s2:wc+s2;
        ry = hc-s2:hc+s2;
        z(rx,ry)=0;    
        %
        [x,y] = find(z==max(z(:)));
        %
        x=x(1);
        y=y(1);
        %
        x_shift = - (x-wc);
        y_shift = - (y-hc);
        
        tform = affine2d([1 0 0; 0 1 0; y_shift x_shift 1]);        
end

