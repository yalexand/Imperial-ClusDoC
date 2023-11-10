function [g_fit, L, N, n, fval] = fit_gr_1c(r,g_exp,N_locs,Area,mode)

                if strcmp(mode,'Gaussian')
                    F = @p_Gauss;
                elseif strcmp(mode,'exponential')
                    F = @p_exp;
                end      
                
    p0 = [30 1000];    
    [x,fval] = run_fitting_gr_1_COMP(r,g_exp,F,Area,p0);
    L = x(1);
    n = x(2);
    N = N_locs/n;

    g_fit = gr_1_COMP(r,F,x,Area);

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
function gr = gr_1_COMP(r,F,p,Area)
    %
    L = p(1);
    n = p(2);
    gr = 1 - 1/n + Area/n*F(r,L);
end
%------------------------------------------------------------------
function [x,fval] = run_fitting_gr_1_COMP(r,g_exp,F,A,p0)
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[5 1 ],[2000 inf],options);
    function y = objfun(p)
        g_theor = gr_1_COMP(r,F,p,A);
        y = norm(g_theor - g_exp);
    end
end
