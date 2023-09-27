function [g_fit, A, L, n, N, F, fval] = fit_gr_14092023(r,g_exp,N_locs,Area,mode,PSF_sigma)

                F = [];
                ro = N_locs/Area;  
                
                fixed_PSF_sigma = ~isempty(PSF_sigma)  && isnumeric(PSF_sigma);                

                if strcmp(mode,'2-comp Gauss')                                   
                    F = {'Gauss','Gauss'};
                elseif strcmp(mode,'2-comp exp')
                    F = {'exp','exp'};                    
                elseif strcmp(mode,'2-comp mixed')
                    F = {'Gauss','exp'};
                % 3-comp    
                elseif strcmp(mode,'Gauss + Gauss Gauss')
                    F = {'Gauss','Gauss','Gauss'};                                        
                elseif strcmp(mode,'Gauss + exp exp')
                    F = {'Gauss','exp','exp'};
                elseif strcmp(mode,'Gauss + mixed')
                    F = {'Gauss','Gauss','exp'};                                        
                end
                %
                if contains(mode,'2-comp')
                    [x,fval,F] = run_fitting_gr_2_COMP(r,g_exp,F);
                    g_fit = gr_2_component(r,F,x);
                    L = x(1:2);
                    A = x(3:4);
                    % estimates presuming equal localizations density                    
                    eta_fit = L(2)/L(1);
                    B_fit = A/Area;
                    denom = B_fit(1)*eta_fit^2 + B_fit(2);
                    n1_fit = B_fit(1)*eta_fit^4/denom^2;
                    n2_fit = B_fit(2)/denom^2;
                    N1_fit = N_locs/(n1_fit + n2_fit*eta_fit^2);
                    N2_fit = N1_fit*eta_fit^2;
                    n = [n1_fit n2_fit];
                    N = [N1_fit N2_fit];
                else
                %
                fixed_PSF_sigma = ~isempty(PSF_sigma)  && isnumeric(PSF_sigma);
                %
                if contains(mode,'+') && ~fixed_PSF_sigma 
                    [x,fval,F] = run_fitting_gr_3_COMP(r,g_exp,F);
                else
                    [x,fval,F] = run_fitting_gr_3_COMP_FPSF(r,g_exp,F,PSF_sigma);
                end
                    g_fit = gr_3_component(r,F,x);                
                    L = x(1:3);
                    A = x(4:6);
                    eta_fit_1 = L(2)/L(1);
                    eta_fit_2 = L(3)/L(1);
                    %
                    B_fit = A/Area;
                    %
                    a = eta_fit_1^2;
                    b = eta_fit_2^2;
                    %
                    denom = a*(b*B_fit(1) + B_fit(3)) + b*B_fit(2);

                    n1_fit = B_fit(1)*a^2*b^2/denom^2;
                    n2_fit = B_fit(2)*b^2/denom^2;
                    n3_fit = B_fit(3)*a^2/denom^2;
                    N1_fit = N_locs/(n1_fit + n2_fit*a + n3_fit*b);
                    N2_fit = N1_fit*a;
                    N3_fit = N1_fit*b;
                    n = [n1_fit n2_fit n3_fit];
                    N = [N1_fit N2_fit N3_fit];
                end
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
function gr = gr_2_component(r,F,p)
    if strcmp(F{1},'Gauss')
        fun1 = @p_Gauss;
    else
        fun1 = @p_exp;
    end
    if strcmp(F{2},'Gauss')
        fun2 = @p_Gauss;
    else
        fun2 = @p_exp;
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
        g_theor = gr_2_component(r,F_in,p);
        y = norm(g_theor - g_exp);
        %y = norm((g_theor - g_exp)./g_exp); %?
    end
    [~,ind]=sort(x(1:2));
    x = x([ind (ind+2)]);
    if ~strcmp(F_in{1},F_in{2})
        F_out = F_in(ind);
    end
end
% 3 component
%------------------------------------------------------------------
function gr = gr_3_component(r,F,p)
    fun1 = @p_Gauss;
    if strcmp(F{2},'Gauss')
        fun2 = @p_Gauss;
    else
        fun2 = @p_exp;
    end
    if strcmp(F{3},'Gauss')
        fun3 = @p_Gauss;
    else
        fun3 = @p_exp;
    end    
    %
    L1 = p(1);
    L2 = p(2);
    L3 = p(3);    
    A1 = p(4);
    A2 = p(5);
    A3 = p(6);
    gr = 1 + A1*fun1(r,L1) + A2*fun2(r,L2) + A3*fun3(r,L3);
end
%------------------------------------------------------------------
function [x,fval,F_out] = run_fitting_gr_3_COMP(r,g_exp,F_in)
    F_out = F_in;
    p0 = [20 50 100 1000 1000 1000];    
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[5 5 5 0 0 0],[50 2000 2000 inf inf inf],options);
    function y = objfun(p)
        g_theor = gr_3_component(r,F_in,p);
        y = norm(g_theor - g_exp);        
    end
    [~,ind]=sort(x(1:3));
    x = x([ind (ind+3)]);
    if ~strcmp(F_in{2},F_in{3})
        F_out = F_in(ind);
    end
end
%------------------------------------------------------------------
function [x,fval,F_out] = run_fitting_gr_3_COMP_FPSF(r,g_exp,F_in,PSF_sigma)
    F_out = F_in;
    p0 = [20 100 1000 1000 1000];    
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[5 5 0 0 0],[50 2000 inf inf inf],options);
    function y = objfun(p)
        g_theor = gr_3_component(r,F_in,[PSF_sigma p]);
        y = norm(g_theor - g_exp);        
    end
    x = [PSF_sigma x];
    [~,ind]=sort(x(1:3));
    x = x([ind (ind+3)]);
    if ~strcmp(F_in{2},F_in{3})
        F_out = F_in(ind);
    end
end

