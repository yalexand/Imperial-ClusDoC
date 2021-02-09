% function [spatialStats, Kr, Lr, Hr, Gr, Xr, tMat] = spatialStats(CroppedPos, Ch1, Ch2, nCores, Conditions)
function [spatialStats, Kr, Lr, Hr, Gr, Xr, tMat] = spatialStats(obj,roi_index,save_dir,Ch1,Ch2,~)

%==========================================================================
%                              FUNCTION
% Function to perform a cross-Ripleys analysis on two overlapping images.
% Based on Ripley's K and derivatives:
%  -Kiskowski et al (2009).  Biophysical Journal 97, pp 1095-1103
%
%==========================================================================
%                               WARNING
% This function is not symmetrical.  K(r) for Im1:Im2 is different than 
% K(r) for Im2:Im1.
%
%==========================================================================
%                             INPUT/OUTPUT
% Inputs:
%   -CroppedPos: Position file produced by LoadCrop.m
%   -Ch1/2: Channels in CroppedIm to analyse; remember, function is not 
%           symmetrical (e.g. Ch1:Ch2 produces a different result than
%           Ch2:Ch1)
%   -nCores: number of cores to use for parallel processing. 1 =
%            non-parallel; 0 = default configuration.
%   -Conditions: Input array containing:
%       -.Bound: maximum distance to perform analysis, in nm
%
%  *** If Conditions is not called in the input, the following input values
%  are used:

DEBUG = false;
if DEBUG
    display ('Using in-file conditions');
    %Ripleys Calculations
    % Conditions.Bound = 1000; %maximum distance to conduct analysis, in nm
    Conditions.Bound = 500; %maximum distance to conduct analysis, in nm
    Conditions.spatialGraph = 1; %graph all spatial data on one graph
    Conditions.RDFgraph = 1; %graph RDF data on its own graph
    Conditions.RipleysGraph = 1; %graph Ripley's K/L/H data on its own graph
    Conditions.fName = 'spatialStats';
else
    Conditions = obj.MIiSR_Conditions;
end

% Outputs:
%   -.Kr/.Lr/.Hr: 2-column array of Ripleys K/L/H function (col 1) and 
%              derivatives (col 2)
%   -.Gr: 1-column array of th Radial Distribution Function.
%   -.Xr: radi (r) values used for calculating K/L/H(r)
%   -.lambda: mean molecular density in the primary image channel
%   -.distanceTable: distance table usd to generate K/H/L/G functions
%
%==========================================================================
%                            DEPENDENCIES
%   Dependencies: A recent copy of Matlab (2010a or newer) plus the 
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

%% Check inputs, convert to table & scale
% commented out - YA
% if (nargin ~= 4) && (nargin ~= 5)
%     edit ('crossRipleys.m');
%     error ('Incorrect number of input arguments.');
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
% cMax = max ([Ch1 Ch2]);
% if cMax > allPos.nCh
%     error ('Ch1 or Ch2 value is greater than the number of available channels');
% end
% 
% % load x/y coords of all molecules
% Im1Ind = allPos.Ch(Ch1).Pos(:,1:3);
% Im2Ind = allPos.Ch(Ch2).Pos(:,1:3);
% 
% clear allPos
% commented out - YA

allPos1 = obj.get_ROI_data_MIiSR(roi_index,Ch1);
allPos2 = obj.get_ROI_data_MIiSR(roi_index,Ch2);
Im1Ind = allPos1(:,1:3);
Im2Ind = allPos2(:,1:3);

nCores = 0;

% start parallel processing
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

%% Pre-processing
% Define bounding box
box(1) = min(min(Im1Ind(:,1)), min(Im2Ind(:,1)));
box(2) = max(max(Im1Ind(:,1)), max(Im2Ind(:,1)));
box(3) = min(min(Im1Ind(:,2)), min(Im2Ind(:,2)));
box(4) = max(max(Im1Ind(:,2)), max(Im2Ind(:,2)));

ImMinSize = min([box(2)-box(1), box(4)-box(3)]);
if ((Conditions.Bound/ImMinSize) > 0.1)
   display ('WARNING: The variable ''Conditions.Bound'' is greater than 10% the size of the image.  This can lead to unusual analysis results.  Consider decreasing the size of Conditions.Bound.');
end

N = size(Im1Ind,1);
lambda = N/((box(2)-box(1))*(box(4)-box(3))); %particle density

%edge correction, particles within Bound units of image's edge
Im1Crop = Im1Ind(Im1Ind(:,1)>(box(1)+Conditions.Bound) & Im1Ind(:,1)<(box(2)-Conditions.Bound) & Im1Ind(:,2)>(box(3)+Conditions.Bound) & Im1Ind(:,2)<(box(4)-Conditions.Bound),:);

%prepare bins & vars
BinW = roundn(Conditions.Bound/100, 1);
Xr = 0:BinW:(Conditions.Bound+BinW); %prepare histograms with one extra bin (to catch stuff > Bound), bins same size as pixels

%used to prepare histogram
tMat = zeros(1,length(Xr));

%% Prepare Distance Tables
parfor ii=1:size(Im1Crop,1)
    pXY = Im1Crop(ii,:)';
    PosList = Im2Ind(Im2Ind(:,1)>=(pXY(1)-Conditions.Bound) & Im2Ind(:,1)<=(pXY(1)+Conditions.Bound) & Im2Ind(:,2)>=(pXY(2)-Conditions.Bound) & Im2Ind(:,2)<=(pXY(2)+Conditions.Bound),1:3)';
    if ~isempty (PosList)
        aa=sum(pXY.*pXY,1); 
        bb=sum(PosList.*PosList,1); 
        ab=pXY'*PosList;
        dTable = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
        tmpHist = hist(dTable,Xr);
        tMat = tMat + tmpHist;
    end
end

tMat = tMat(:,1:(end-1)); %remove last bin

for kk=1:length(tMat)
    Kr(kk,1) = sum(tMat(1:kk))/N;
end

%K, L and H Calculation
Kr = Kr./lambda;
Lr = sqrt(Kr./pi());
for ii=1:length(Xr)-1
    Hr(ii,1)=Lr(ii)-(Xr(ii));  
end

%K', L' and H' Calculation (derivatives of K,L & H)
for ii=2:(length(Xr)-2)
	Kr(ii-1,2) = (Kr(ii+1)-Kr(ii-1))/(BinW*2);
    Lr(ii-1,2) = (Lr(ii+1)-Lr(ii-1))/(BinW*2);
    Hr(ii-1,2) = (Hr(ii+1)-Hr(ii-1))/(BinW*2);
end

%G(r) (RDF) Calculation
for ii=1:length(tMat)
    Gr(ii) = tMat(ii)/(2*pi()*lambda*BinW*(Xr(ii+1))); %2D
end

%export data to structure
spatialStats.Kr = Kr;
spatialStats.Lr = Lr;
spatialStats.Hr = Hr;
spatialStats.Gr = Gr';
spatialStats.Xr = (Xr(1:end-1))';
spatialStats.lambda = lambda;
spatialStats.distanceTable = tMat;

mkdir(save_dir,['roi_' num2str(roi_index)]);
roi_save_dir = [save_dir filesep 'roi_' num2str(roi_index)];

%% Graph Data
h = figure(1);
set (figure(1), 'Color', 'w');
subplot(2,2,1);
kMax = ceil(max(Kr(:,1))*1.1);
plot(spatialStats.Xr, Kr(:,1),'k');
axis([0 Conditions.Bound, 0, kMax]);
xlabel ('Radius (nm)');
ylabel ('K(r)');
title ('Ripley''s K-function');

subplot(2,2,2);
kMax = ceil(max(Lr(:,1))*1.1);
plot(spatialStats.Xr, Lr(:,1),'k');
axis([0 Conditions.Bound, 0, kMax]);
xlabel ('Radius (nm)');
ylabel ('L(r)');
title ('Ripley''s L-function');

subplot(2,2,3);
kMax = ceil(max(Hr(:,1))*1.1);
kMin = floor(min(Hr(:,1))*1.1);
plot(spatialStats.Xr, Hr(:,1),'k');
axis([0 Conditions.Bound, kMin, kMax]);
xlabel ('Radius (nm)');
ylabel ('H(r)');
title ('Ripley''s H-function');

subplot(2,2,4);
kMax = ceil(max(Gr)*1.1);
kMin = floor(min(Gr)*1.1);
plot(spatialStats.Xr, Gr','k');
axis([0 Conditions.Bound, kMin, kMax]);
xlabel ('Radius (nm)');
ylabel ('RDF/G(r)');
title ('Radial Distribution Function');

tPos = get(figure(1), 'Position');
set(figure(1), 'Position', [(tPos(1)/2), (tPos(2)/2), 1120, 840]);
saveas (figure(1), [roi_save_dir filesep 'Spatial Stats - ' Conditions.fName '.fig']);
% saveas (figure(1), ['Spatial Stats - ' Conditions.fName '.pdf']);
pause(0.5);
close (h);

h = figure(2);
set (figure(2), 'Color', 'w');
kMax = ceil(max(Kr(:,1))*1.1);
plot(spatialStats.Xr, Kr(:,1),'k');
axis([0 Conditions.Bound, 0, kMax]);
xlabel ('Radius (nm)');
ylabel ('K(r)');
title ('Ripley''s K-function');
%saveas (figure(2), ['K-Function - ' Conditions.fName '.fig']);
pause(0.5);
close (h);

h = figure(3);
kMax = ceil(max(Lr(:,1))*1.1);
plot(spatialStats.Xr, Lr(:,1),'k');
axis([0 Conditions.Bound, 0, kMax]);
xlabel ('Radius (nm)');
ylabel ('L(r)');
title ('Ripley''s L-function');
% saveas (figure(3), ['L-Function - ' Conditions.fName '.fig']);
pause(0.5);
close (h);

h = figure(4);
kMax = ceil(max(Hr(:,1))*1.1);
kMin = floor(min(Hr(:,1))*1.1);
plot(spatialStats.Xr, Hr(:,1),'k');
axis([0 Conditions.Bound, kMin, kMax]);
xlabel ('Radius (nm)');
ylabel ('H(r)');
title ('Ripley''s H-function');
saveas (figure(4), [roi_save_dir filesep 'H-Function - ' Conditions.fName '.fig']);
pause(0.5);
close (h);

h = figure(5);
kMax = ceil(max(Gr)*1.1);
kMin = floor(min(Gr)*1.1);
plot(spatialStats.Xr, Gr','k');
axis([0 Conditions.Bound, kMin, kMax]);
xlabel ('Radius (nm)');
ylabel ('RDF/G(r)');
title ('Radial Distribution Function');
saveas (figure(5), [roi_save_dir filesep 'RDF - ' Conditions.fName '.fig']);
pause(0.5);
close (h);

%% Clean-up
if openPool == 1 %close pool of opened by function
    delete (gcp);
end

end %end function

function outnum = roundn(innum,varargin)
% Round numbers to a specified power of 10.
% Syntax
% outnum = roundn(innum)
% outnum = roundn(innum,n)
% Description
% outnum = roundn(innum) rounds the elements of innum to the nearest one-hundredth.
% outnum = roundn(innum,n) specifies the power of 10 to which the elements of innum are rounded. For example, if n = 2, round to the nearest hundred (10^2). 

numvarargs = length(varargin);
if numvarargs > 1
    error('roundn','requires at most 1 optional input');
end

optargs = {-2};

% skip any new inputs if they are empty
newVals = ~cellfun('isempty', varargin);
optargs(newVals) = varargin(newVals);

% Place optional args in named variables
[n] = optargs{:};
factor = 10^(-n);
outnum = round(innum*factor)/factor;
end