function [g_fit, A, L, n, p, N, F, fval] = fit_gr3(r,g_exp,N_locs,Area,mode,jointA,jointAA,PSF_sigma)

A1 = []; 
A2 = []; 
A3 = []; 
L1 = []; 
L2 = []; 
L3 = []; 
p1 = []; 
p2 = []; 
p3 = []; 
                ro = N_locs/Area;  
                
                fixed_PSF_sigma = ~isempty(PSF_sigma)  && isnumeric(PSF_sigma);                

                F1 = @p_Gauss;
                %
                if strcmp(mode,'Gaussian') 
                    F2 = @p_Gauss;
                    F3 = @p_Gauss;
                    F = {'Gauss','Gauss','Gauss'};
                elseif strcmp(mode,'exponential')
                    F2 = @p_exp;
                    F3 = @p_exp;
                    F = {'Gauss','exp','exp'};
                elseif strcmp(mode,'mixed')
                    F2 = @p_Gauss;
                    F3 = @p_exp;
                    F = {'Gauss','Gauss','exp'};
                end      

                if ~jointAA && ~fixed_PSF_sigma && ~jointA          % 8 - remove ?
                    x0 = [1 1 1 7 30 300 0.3 0.3];
                    [x,fval] = run_fitting_gr_NOT_FPCFD_NOT_JA(r,g_exp,F1,F2,F3,x0);
                        A1 = x(1);
                        A2 = x(2);
                        A3 = x(3);                        
                        L1 = x(4);
                        L2 = x(5);
                        L3 = x(6);                        
                        p1 = x(7);
                        p2 = x(8);                                      
                elseif ~jointAA && ~fixed_PSF_sigma && jointA
                    x0 = [1 1 7 30 300 0.3 0.3];
                    [x,fval] = run_fitting_gr_NOT_FPCFD_AND_JA(r,g_exp,F1,F2,F3,x0);
                        A1 = x(1);
                        A2 = x(2);
                        A3 = x(2);                        
                        L1 = x(3);
                        L2 = x(4);
                        L3 = x(5);                        
                        p1 = x(6);
                        p2 = x(7);                        
                elseif ~jointAA && fixed_PSF_sigma && ~jointA           % 7 - remove ?
                    x0 = [1 1 1 7 300 0.3 0.3];
                    [x,fval] = run_fitting_gr_AND_FPCFD_NOT_JA(r,g_exp,F1,F2,F3,PSF_sigma,x0);
                        A1 = x(1);
                        A2 = x(2);
                        A3 = x(3);                        
                        L1 = PSF_sigma;
                        L2 = x(4);
                        L3 = x(5);                        
                        p1 = x(6);
                        p2 = x(7);                        
                elseif ~jointAA && fixed_PSF_sigma && jointA
                    x0 = [1 1 7 300 0.3 0.3];
                    [x,fval] = run_fitting_gr_AND_FPCFD_AND_JA(r,g_exp,F1,F2,F3,PSF_sigma,x0);
                        A1 = x(1);
                        A2 = x(2);
                        A3 = x(2);                        
                        L1 = PSF_sigma;
                        L2 = x(3);
                        L3 = x(4);                        
                        p1 = x(5);
                        p2 = x(6);                        
                elseif jointAA && fixed_PSF_sigma 
                    x0 = [1 1 1 30 300];
                    [x,fval] = run_fitting_gr_AND_FPCFD_AND_JA_ALL(r,g_exp,F1,F2,F3,PSF_sigma,x0);
                        A1 = x(1);
                        A2 = x(2);
                        A3 = x(3);                        
                        L1 = PSF_sigma;
                        L2 = x(4);
                        L3 = x(5);                        
                        A = A1 + A2 + A3;
                        p1 = A1/A;
                        p2 = A2/A;
                        A1 = A;
                        A2 = A;
                        A3 = A; 
                elseif jointAA && ~fixed_PSF_sigma 
                    x0 = [1 1 1 7 30 300];
                    [x,fval] = run_fitting_gr_NOT_FPCFD_AND_JA_ALL(r,g_exp,F1,F2,F3,x0);
                        A1 = x(1);
                        A2 = x(2);
                        A3 = x(3);                        
                        L1 = x(4);
                        L2 = x(5);
                        L3 = x(6);                        
                        A = A1 + A2 + A3;
                        p1 = A1/A;
                        p2 = A2/A;
                        A1 = A;
                        A2 = A;
                        A3 = A; 
                end
                
                p3 = 1 - p1 - p2;
                
                g_fit = p1*(1 + A1*F1(r,L1)) + p2*(1 + A2*F2(r,L2)) + p3*(1 + A3*F3(r,L3));
                                
L = [L1 L2 L3];
A = [A1 A2 A3];
p = [p1 p2 p3];
[~,sortind] = sort(L);
%
L = L(sortind);
A = A(sortind);
p = p(sortind);
F = F(sortind);
%
n = 2*pi.*A*ro.*L.^2;
N = N_locs*p./n;

%                 disp('start!')                
%                 L
%                 p
%                 n
%                 F
%                 disp('bingo!')  

end
%------------------------------------------------------------------
function f = p_exp(r,L)
    f = exp(-r/L);
end
%------------------------------------------------------------------
function f = p_Gauss(r,L)
    f = 1/2*exp(-r.^2/(4*L^2));
end

%------------------------------------------------------------------
function [x,fval] = run_fitting_gr_NOT_FPCFD_NOT_JA(r,g_exp,fun1,fun2,fun3,p0)
% convention
%     A1
%     A2
%     A3
%     L1
%     L2
%     L3
%     p1
%     p2
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[0 0 0 5 5 5 0 0],[inf inf inf 2000 2000 2000 1 1],options);
    function y = objfun(p)
        g1 = 1 + p(1)*fun1(r,p(4));
        g2 = 1 + p(2)*fun2(r,p(5));
        g3 = 1 + p(3)*fun3(r,p(6));        
        g_theor = p(7)*g1 + p(8)*g2 + (1-p(7)-p(8))*g3;
        y = norm(g_theor - g_exp);
    end
end

%------------------------------------------------------------------
function [x,fval] = run_fitting_gr_NOT_FPCFD_AND_JA(r,g_exp,fun1,fun2,fun3,p0)
% convention
%     A1
%     A
%     L1
%     L2
%     L3
%     p1
%     p2
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[0 0 5 5 5 0 0],[inf inf 2000 2000 2000 1 1],options);
    function y = objfun(p)
        g1 = 1 + p(1)*fun1(r,p(3));
        g2 = 1 + p(2)*fun2(r,p(4));
        g3 = 1 + p(2)*fun3(r,p(5));        
        g_theor = p(6)*g1 + p(7)*g2 + (1-p(6)-p(7))*g3;
        y = norm(g_theor - g_exp);
    end
end

%------------------------------------------------------------------
function [x,fval] = run_fitting_gr_AND_FPCFD_NOT_JA(r,g_exp,fun1,fun2,fun3,PSF_sigma,p0)
% convention
%     A1
%     A2
%     A3
%     missing - fixed
%     L2
%     L3
%     p1
%     p2
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[0 0 0 5 5 0 0],[inf inf inf 2000 2000 1 1],options);
    function y = objfun(p)
        g1 = 1 + p(1)*fun1(r,PSF_sigma);
        g2 = 1 + p(2)*fun2(r,p(4));
        g3 = 1 + p(3)*fun3(r,p(5));        
        g_theor = p(6)*g1 + p(7)*g2 + (1-p(6)-p(7))*g3;
        y = norm(g_theor - g_exp);
    end
end

%------------------------------------------------------------------
function [x,fval] = run_fitting_gr_AND_FPCFD_AND_JA(r,g_exp,fun1,fun2,fun3,PSF_sigma,p0)
% convention
%     A1
%     A
%     missing - fixed
%     L2
%     L3
%     p1
%     p2
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[0 0 5 5 0 0],[inf inf 2000 2000 1 1],options);
    function y = objfun(p)
        g1 = 1 + p(1)*fun1(r,PSF_sigma);
        g2 = 1 + p(2)*fun2(r,p(3));
        g3 = 1 + p(2)*fun3(r,p(4));        
        g_theor = p(5)*g1 + p(6)*g2 + (1-p(5)-p(6))*g3;
        y = norm(g_theor - g_exp);
    end
end

%------------------------------------------------------------------
function [x,fval] = run_fitting_gr_NOT_FPCFD_AND_JA_ALL(r,g_exp,fun1,fun2,fun3,p0)
% convention
%     A1
%     A2
%     A3
%     L1
%     L2
%     L3
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[0 0 0 5 5 5],[inf inf inf 2000 2000 2000],options);
    function y = objfun(p)
        A1 = p(1);
        A2 = p(2);
        A3 = p(3);
        A = A1 + A2 + A3;
        p1 = A1/A;
        p2 = A2/A;
        p3 = A3/A;
        g1 = 1 + A*fun1(r,p(4));
        g2 = 1 + A*fun2(r,p(5));
        g3 = 1 + A*fun3(r,p(6));  
        %                
        g_theor = p1*g1 + p2*g2 + p3*g3;
        y = norm(g_theor - g_exp);
    end
end

%------------------------------------------------------------------
function [x,fval] = run_fitting_gr_AND_FPCFD_AND_JA_ALL(r,g_exp,fun1,fun2,fun3,PSF_sigma,p0)
% convention
%     A1
%     A2
%     A3
%     fixed
%     L2
%     L3
    Nmax = 10000*length(p0);
    options = optimset('MaxFunEvals',Nmax,'MaxIter',Nmax);
    [x,fval] = fminsearchbnd(@objfun,p0,[0 0 0 5 5],[inf inf inf 2000 2000],options);
    function y = objfun(p)
        A1 = p(1);
        A2 = p(2);
        A3 = p(3);
        A = A1 + A2 + A3;
        p1 = A1/A;
        p2 = A2/A;
        p3 = A3/A;
        g1 = 1 + A*fun1(r,PSF_sigma);
        g2 = 1 + A*fun2(r,p(4));
        g3 = 1 + A*fun3(r,p(5));  
        %                
        g_theor = p1*g1 + p2*g2 + p3*g3;
        y = norm(g_theor - g_exp);
    end
end
