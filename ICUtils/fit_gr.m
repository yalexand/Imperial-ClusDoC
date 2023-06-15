function [g_fit, L1 L2 n1 n2 p1 p2 N1 N2, fval] = fit_gr(r,g_exp,N_locs,Area,mode,num_param)

                ro = N_locs/Area;  
                
                if strcmp(mode,'Gaussian')
                    F1 = @p_Gauss;
                    F2 = @p_Gauss;
                elseif strcmp(mode,'exponential')
                    F1 = @p_exp;
                    F2 = @p_exp;
                elseif strcmp(mode,'mixed')
                    F1 = @p_Gauss;
                    F2 = @p_exp;
                end      
                
                if 5 == num_param
                        [x,fval] = run_fitting_gr_5p(r,g_exp,F1,F2,[1 30 1 100, 0.5]);
                        A1 = x(1);
                        L1 = x(2);
                        A2 = x(3);
                        L2 = x(4);
                        p1 = x(5);
                        p2 = 1 - p1;
                elseif 4 == num_param
                        [x,fval] = run_fitting_gr_4p(r,g_exp,F1,F2,[1 30 100, 0.5]);
                        A = x(1);
                        L1 = x(2);
                        L2 = x(3);
                        p1 = x(4);
                        p2 = 1 - p1;
                        A1 = A;
                        A2 = A;
                end                   
                        g_fit = p1*(1 + A1*F1(r,L1)) + p2*(1 + A2*F2(r,L2));
                     
                        L_min = min(L1,L2);
                        L_max = max(L1,L2);
                        if L_min==L1
                            A_min = A1;
                            A_max = A2;
                            p_min = p1;
                        else
                            A_min = A2;
                            A_max = A1;
                            p_min = p2;
                        end

                        A1 = A_min;
                        L1 = L_min;
                        A2 = A_max;
                        L2 = L_max;
                        p1 = p_min;
                        p2 = 1 - p_min;

                        n1 = 2*pi*A1*ro*L1^2;
                        n2 = 2*pi*A2*ro*L2^2;

                        N1 = N_locs*p1/n1;
                        N2 = N_locs*p2/n2;
                        %N_locs - (N1*n1 + N2*n2)                         
end
%------------------------------------------------------------------
function f = p_exp(r,L)
    f = exp(-r/L);
end
%------------------------------------------------------------------
function f = p_Gauss(r,L)
    f = sqrt(2/pi)*exp(-r.^2/(2*L^2));
end
%------------------------------------------------------------------
function [x,fval] = run_fitting_gr_4p(r,g_exp,fun1,fun2,p0)
% convention
%     A
%     L1
%     L2
%     p1
    Nmax = 400*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[0 3 3 0],[inf 1000 1000 1],options);
    function y = objfun(p)
        g1 = 1 + p(1)*fun1(r,p(2));
        g2 = 1 + p(1)*fun2(r,p(3));
        g_theor = p(4)*g1 + (1-p(4))*g2;
        y = norm(g_theor - g_exp);
    end
end
%------------------------------------------------------------------
function [x,fval] = run_fitting_gr_5p(r,g_exp,fun1,fun2,p0)
% convention
%     A1
%     L1
%     A2
%     L2
%     p1
    Nmax = 400*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[0 3 0 3 0],[inf 1000 inf 1000 1],options);
    function y = objfun(p)
        g1 = 1 + p(1)*fun1(r,p(2));
        g2 = 1 + p(3)*fun2(r,p(4));
        g_theor = p(5)*g1 + (1-p(5))*g2;
        y = norm(g_theor - g_exp);
    end
end
