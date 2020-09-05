%-------------------------------------
function u_out = render_to_nbits(nbits,u,qthresh)
    s = u(u>0);
    t = quantile(s(:),qthresh);
    u(u>t)=t;
    u_out = map(u,0,2^nbits-1);    
end
