function [OPTICStable, RD, CD, order, Epsilon] = OPTICSStats(obj,roi_index,chan,save_dir,~)

%==========================================================================
%                              FUNCTION
% Performs OPTICS analysis for cluster identification on molecule position
% file.
%==========================================================================
%                             INPUT/OUTPUT
% Inputs:
%   -CroppedPos: Position file produced by GSDLoadCrop
%   -Ch1: Which channel in CroppedPos to analyse
%   -nCores: number of cores to use for parallel processing. 1 =
%            non-parallel; 0 = default configuration.
%   -k = minimum number of particles to be considered a cluster.
%
% Outputs:
%   -RD = reachability distances for all particles (m,1)
%   -CD = table of core distances
%   -order = vector providing the order of the particles in n.
%   -OPTICStable: table of:
%       -Col 1: order
%       -Col 2: RD
%       -Col 3: CD
%       -Col 4/5/6: X/Y/Z coordinates
%   -Epsilon = Reachability distance of randomly distributed data
%
%==========================================================================
%                            DEPENDENCIES
%   Dependencies: A recent copy of Matlab (2010a or newer) plus the 
%                 parallel processing toolbox. This code is processor-
%                 intensive, thus a high-end, multi-processor workstation 
%                 is recommended.
%
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

% %% Check Inputs
% if nargin ~= 4
%     edit ('OPTICS.m');
%     error ('OPTICS requires 4 input variables.');
% end
% 
% if k < 2
%     display ('Minimum cluster size is 2.  K is now set as 2');
%     k = 2;
% end
% 
% if ~(exist (CroppedPos, 'file'))
%     error ('Input file is not in working directory');
% end
% 
% load (CroppedPos);
% if exist ('cropPos', 'var')
%     allPos = cropPos;
%     clear cropPos
% end
% 
% if Ch1 > allPos.nCh
%     error ('Ch1 value is greater than the number of available channels');
% end
% 
% %Extract x/y/z
% pXY = allPos.Ch(Ch1).Pos(:,1:3);
% clear allPos;

if obj.MIiSR_Specs.densityFilter(chan)
    allPos = obj.densityFilter(roi_index,chan);
else
    allPos = obj.get_ROI_data_MIiSR(roi_index,chan);
end

pXY = allPos(:,1:3);
clear('allPos');

nCores = 0;

k = obj.MIiSR_Specs.OPTICSk(chan);

m = size(pXY,1);

if max(pXY(:,3)) == 0 %if 2D array
    Epsilon = ((prod(max(pXY(:,1:2))-min(pXY(:,1:2)))*k*gamma(.5*1+1))/(m*sqrt(pi.^2))).^(1/2);
else
    Epsilon = ((prod(max(pXY)-min(pXY))*k*gamma(.5*1+1))/(m*sqrt(pi.^2))).^(1/2);
end

%Prepare for parallel processing
openPool = 0;
poolInfo = gcp('nocreate');

if isempty(poolInfo) && (nCores ~= 1) %if no pool is open & parallel processing is selected
    openPool = 1;
    if (nCores == 0)
        parpool; %open default settings
    else
        parpool (nCores); %open set number of cores
    end
end

tic

%%  Generate CD, RD & order tables
RD = zeros(m,1);
CD = RD;

parfor ii=1:m
    tPos = pXY(:,1:3); %extract coords
    tPos(ii,:) = []; %eliminate self-position
    D = sort(sqrt((tPos(:,1)-pXY(ii,1)).^2 + (tPos(:,2)-pXY(ii,2)).^2 + (tPos(:,3)-pXY(ii,3)).^2));
    CD(ii)=D(k-1); %CD as most distant particle in local cluster of size k
end

% Calculate Core Distances
pXY(:,4) = 1:1:m; %col 4 becomes order ID number for particles
pXY(:,5) = CD;

clear CD

%calculate RD
for ii=1:(m-1)
    [RD(ii+1), nnD] = min(sqrt((pXY(ii+1:end,1)-pXY(ii,1)).^2 + (pXY(ii+1:end,2)-pXY(ii,2)).^2 + (pXY(ii+1:end,3)-pXY(ii,3)).^2)); %distance to all subsiquent points
    tPos = pXY(ii+1,:); %prepare to rearrange dataset
    pXY(ii+1,:) = pXY(ii+nnD,:); %move closest particle to below current particle
    pXY(ii+nnD,:) = tPos;
    RD(ii+1) = max([RD(ii+1), pXY(ii+1,5)]); %set RD
end

%calculate RD for particle 1; set as Epsilon *2
RD(1)=Epsilon*2;
CD = pXY(:,5);

%write output
order = pXY(:,4);
OPTICStable(:,1) = order;
OPTICStable(:,2) = RD;
OPTICStable(:,3) = CD;
OPTICStable(:,4:6) = pXY(:,1:3);

%% Clean-up
% if openPool == 1 %close pool of opened by function
%     delete (gcp);
% end

disp(['OPTICS roi ' num2str(roi_index) ', channel ' num2str(chan) ',time ' num2str(toc)]);

end
