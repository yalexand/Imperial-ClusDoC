function sigma = SMLM_sigma_calculator(NA,lambda,pix_size)

    sigma = Rayleigh_resolution(lambda,NA)/pix_size/2;
    
end

%
function r = Rayleigh_resolution(lambda,NA)
    r = 0.61*lambda/NA;
end