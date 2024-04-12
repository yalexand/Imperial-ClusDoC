function [g_fit, L1, L2, N1, N2, n1, n2, F_out, fval] = fit_gr_2c_5p(r,g_exp,N_locs,Area,mode)

    F = strsplit(mode,' ');
                
    rho_avr = N_locs/Area;
    rho_max = 3/16;

    Z_min = 5;
    Z_max = 2000;
    N_min = locs_per_cluster(Z_min,rho_avr);
    N_max = N_locs/2; 
    n_min = 1; 
    n_max = N_locs/n_min;
    
    [x,fval,F] = run_fitting_gr_2_COMP(r,g_exp,F);
    L1 = x(1);
    L2 = x(2);
    
    p0 = [ 30 200 1000 ];
    MIN = [N_min N_min n_min];
    MAX = [N_max N_max n_max];    
    [x,fval,F_out,n12_flipped] = run_fitting_assisted_L1L2(r,g_exp,F,N_locs,Area,p0,MIN,MAX,L1,L2);    
    
                N1 = x(1);
                N2 = x(2);
                if ~n12_flipped
                    n1 = x(3);
                    n2 = (N_locs - n1*N1)/N2; 
                else
                    n2 = x(3);
                    n1 = (N_locs - n2*N2)/N1; 
                            L1
                            L2
                            N1
                            N2
                            n1
                            n2
                end

    g_fit = multi_component_gr(r,F,[L1 L2],x(1:2),x(3),N_locs,Area) ;
end
%------------------------------------------------------------------
function f = p_exp(r,L)
    f = 1/(2*pi*L^2)*exp(-r/L);
end
%------------------------------------------------------------------
function f = p_Gauss(r,L)
    f = 1/(4*pi*L^2)*exp(-r.^2/(4*L^2));
end
%------------------------------------------------------------------
function f = p_disk(r,R) 
    f = 1./(2*pi*r).*( 4*r/(pi*R^2).*acos(r/(2*R)) - 2*r.^2/(pi*R^3).*(1 - r.^2/(4*R^2)).^(1/2) );
    f(imag(f)~=0) = 0;
end
%------------------------------------------------------------------
function gr = gr_2_COMP(r,F,p)
    if strcmp(F{1},'Gauss')
        fun1 = @p_Gauss;        
    elseif strcmp(F{1},'exp')
        fun1 = @p_exp;
    elseif strcmp(F{1},'disk')
        fun1 = @p_disk;        
    end
    if strcmp(F{2},'Gauss')
        fun2 = @p_Gauss;        
    elseif strcmp(F{2},'exp')
        fun2 = @p_exp;
    elseif strcmp(F{2},'disk')
        fun2 = @p_disk;        
    end    
    %
    L1 = p(1);
    L2 = p(2);
    A1 = p(3);
    A2 = p(4);
    gr = 1 + A1*fun1(r,L1) + A2*fun2(r,L2);
end
%------------------------------------------------------------------
function [x,fval,F_out] = run_fitting_gr_2_COMP(r,g_exp,F_in)
    F_out = F_in;
    p0 = [30 100 1000 1000];    
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[5 5 0 0],[2000 2000 inf inf],options);
    function y = objfun(p)
        g_theor = gr_2_COMP(r,F_in,p);
        y = norm(g_theor - g_exp);
        %y = norm((g_theor - g_exp)./g_exp); %?
    end
    [~,ind]=sort(x(1:2));
    x = x([ind (ind+2)]);
    if ~strcmp(F_in{1},F_in{2})
        F_out = F_in(ind);
    end
end
%------------------------------------------------------------------
function gr_theor = multi_component_gr(r,F,Z,N,n_,Ntot,Area) 

    if strcmp(F{1},'Gauss')
        fun1 = @p_Gauss;        
    elseif strcmp(F{1},'exp')
        fun1 = @p_exp;
    elseif strcmp(F{1},'disk')
        fun1 = @p_disk;        
    end
    if strcmp(F{2},'Gauss')
        fun2 = @p_Gauss;        
    elseif strcmp(F{2},'exp')
        fun2 = @p_exp;
    elseif strcmp(F{2},'disk')
        fun2 = @p_disk;        
    end    

    rho = Ntot/Area;
    
    rank = length(N);
    n_rank = max(0,(Ntot - sum(N(1:(rank-1)).*n_))/N(rank));
    n = [n_ n_rank];

    acc = zeros(size(r));
    for c=1:numel(n)
        if 1==c
            fun=fun1;
        else
            fun=fun2;
        end
        acc = acc + n(c)*N(c)^2*fun(r,Z(c));
    end
    gr_theor = 1 + acc*(1/(Area*rho^2));

end
%------------------------------------------------------------------
function [x,fval,F_out, n12_flipped] = run_fitting_assisted_L1L2(r,g_exp,F_in,Ntot,Area,p0,MIN,MAX,L1,L2)
    % Z
    % N
    % n
    F_out = F_in;
    %
    n12_flipped = false;
    %
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,MIN,MAX,options);
    %
    function y = objfun(p)
        g_theor = multi_component_gr(r,F_in,[L1 L2],p(1:2),p(3),Ntot,Area);
        y = norm(g_theor - g_exp);
    end
    [~,ind]=sort(x(1:2));
    x = x([ind 3]);
    %
    if ind(1)==2, n12_flipped = true; end 
    %
    if ~strcmp(F_in{1},F_in{2})
        F_out = F_in(ind);
    end    
end
%------------------------------------------------------------------
function N = locs_per_cluster(Z,rho)
    N = pi*(1.5*Z)^2*rho;
end
