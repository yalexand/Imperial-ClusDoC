function [beta,dbeta] = RDFquant (MIiSRfile, par0)
%==========================================================================
%                              FUNCTION
% Fit's G(r) data to equation 10 of Sengupta et al, Nature Methods 8,
% 969 (2011). Provides an approximation of the number of molecules/cluster
% from RDF plots
%==========================================================================
%                               WARNING
% The function of Sengupta et al does not provide an exact count of
% molecules in a cluster in most super-resolution images. Rather, it counts 
% the total number of fluorophore detections in a sample. As such,
% over-sampling will increase the number of molecules per cluster, and
% under-sampling will reduce the number of molecules per cluster.
%==========================================================================
%                             INPUT/OUTPUT
% -MIiSRfile: Data file (.mat) produced following MIiSR processing. Both 
%              SAA and RDF analysis must have been performed on the data
%              set
% -par0: A 1*2 array containing initial guess for 
%        [peak, molecules/cluster]; default is [1 1]
%
% -beta = [cluster radius, cluster size]
% -dbeta = uncertainties from beta
%==========================================================================
%                            DEPENDENCIES
%   Dependencies: A recent copy of Mtlab (2010a or newer) plus the 
%                 parallel processing toolbox. This code is processor-
%                 intensive, thus a high-end, multi-processor workstation 
%                 is recommended.
%==========================================================================
%                             CITATION
%
% This script is provided as a supplemental material in:
%
%   Fabiana A. Caetano, Brennan S. Dirk, Joshua H.K. Tam, P. Craig 
%       Cavanagh, Maria Goiko, Stephen S.G. Ferguson, Stephen H. Pasternak,
%       Jimmy D. Dikeakos, John R. de Bruyn, Bryan Heit. MIiSR: Analysis of 
%		Molecular Interactions in Super-Resolution Imaging Enables the Study
%		of Protein Interactions, Dynamics and Formation of Multi-protein 
%		Structures. 2015. PLoS Computational Biology
% 
% Please reference this paper in any publications which use this script for 
% analysis.
%==========================================================================

%% Check inputs, set global variables
if nargin == 1
    par0 = [1, 1];
elseif nargin ~= 2
    error ('Incorrect number of input varibales has been called');
end

if ~(exist (MIiSRfile, 'file'))
    error ('Input file is not in working directory');
end

load (MIiSRfile, 'MIiSRdata');

if ~isfield (MIiSRdata, 'SAAdata') && ~isfield (MIiSRdata, 'spatialData')
    error ('Both SAA and RDF anlaysis must be performed on the data set');
end

global sig_s rhoavg % set global variables

%determine precission of primary channel
if MIiSRdata.queueInfo.spatialCh1 == 1
    sig_s = MIiSRdata.SAAdata.Im1MeanPrecission; % mean precission error of primary channel
elseif MIiSRdata.queueInfo.spatialCh1 == 2
    sig_s = MIiSRdata.SAAdata.Im2MeanPrecission; % mean precission error of primary channel
else MIiSRdata.queueInfo.spatialCh1
    sig_s = MIiSRdata.SAAdata.Im3MeanPrecission; % mean precission error of primary channel
end

rhoavg = MIiSRdata.spatialData.lambda; %average density

%fit to equation 10 from paper
[beta,res,j,cov,mse] = nlinfit(MIiSRdata.spatialData.Xr,MIiSRdata.spatialData.Gr,@eq10,par0);

%calculate precission of fit
ci = nlparci(beta,res,'covar',cov);
ci2 = ci';
dbeta = abs(beta - ci2(1,:))/2;

clear sig_s rhoave

end % end main function

function yfit = eq10(p,x)
%
% function for g(r) defined in Eq. (10) of the above paper

% need to enter values for these parameters:
%      sigma_s, width of pt spread function
%          set as global, above;  sig_s = 1;
%      average density of proteins
%          set as global, above;  rhoavg = 1;

global sig_s rhoavg

gpsf = 1./(4*pi*sig_s^2) * exp(-x.^2/(4*sig_s^2));

yfit = (1/rhoavg + p(1) * exp(-x/p(2)) + 1) .* gpsf;

end %end eq10