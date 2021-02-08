function [RDFmap, RDFtable] = RDFsegment (CroppedPos, Ch1, nCores, cRadius)

%%==========================================================================
%                              FUNCTION
% This function uses Ripley's H function at a set radius to produce a 
% density map which can be used for further quantification of clustering.
% This is based on the approach of Owen D.M, Methods in Enzymology, 2012
%
%==========================================================================
%                             INPUT/OUTPUT
% Inputs:
%   -CroppedPos: Position file produced by LoadCrop.m
%   -Ch1: Channel to analyze in CroppedPos
%   -nCores: number of cores to use for parallel processing. 1 =
%            non-parallel; 0 = default configuration.
%   -cRadius: Radius to use when calculating local density, in nm
% Outputs:
%   -RDFmap: 8-bit image, with 20nm pixels, of image density. Is scaled 0 -
%            255
%   -RDFtable: 4-column table containing the H(cRadius) of each molecule in
%              the dataset, as well as the molecules X/Y/Z coordinates.
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
if nargin ~= 4
    error ('Incorrect number of input arguments.');
end

if ~(exist (CroppedPos, 'file'))
    error ('Input file is not in working directory');
end

load (CroppedPos);
if exist ('cropPos', 'var')
    allPos = cropPos;
    clear cropPos
end

if Ch1 > allPos.nCh
    error ('Value entered for Ch1 is greater than the number of available channels');
end

% load x/y/z coords of all molecules
Im1Ind = allPos.Ch(Ch1).Pos(:,1:3);

clear allPos

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
box(1) = min(Im1Ind(:,1));
box(2) = max(Im1Ind(:,1));
box(3) = min(Im1Ind(:,2));
box(4) = max(Im1Ind(:,2));

ImMinSize = min([box(2)-box(1), box(4)-box(3)]);
if ((cRadius/ImMinSize) > 0.1)
   display ('WARNING: The variable ''cRadius'' is greater than 10% the size of the image.  This can lead to unusual analysis results.  Consider decreasing the size of cRadius.');
end

%edge correction, particles within Bound units of image's edge
Im1Crop = Im1Ind(Im1Ind(:,1)>(box(1)+cRadius) & Im1Ind(:,1)<(box(2)-cRadius) & Im1Ind(:,2)>(box(3)+cRadius) & Im1Ind(:,2)<(box(4)-cRadius),:);

%prepare bins & vars
BinW = ceil(cRadius/50); %round up to nearest intiger
Xr = 0:BinW:(cRadius+BinW); %prepare histograms with one extra bin (to catch stuff > Bound), bins same size as pixels

%% Prepare Distance Tables
RDFtable = zeros(size(Im1Crop,1),4);

parfor ii=1:size(Im1Crop,1)
    pXY = Im1Crop(ii,:)';
    PosList = Im1Ind(Im1Ind(:,1)>=(pXY(1)-cRadius) & Im1Ind(:,1)<=(pXY(1)+cRadius) & Im1Ind(:,2)>=(pXY(2)-cRadius) & Im1Ind(:,2)<=(pXY(2)+cRadius),1:3)';
    if ~isempty (PosList)
        aa=sum(pXY.*pXY,1); 
        bb=sum(PosList.*PosList,1); 
        ab=pXY'*PosList;
        dTable = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
        tmpHist = hist(dTable,Xr);
        lambda = sum(tmpHist)/cRadius^2; %loal density
        tmpHist(end) = [];
        
        %calculate single molicule K/L/R values
        Kr = zeros(1,length(tmpHist));
        for kk=1:length(tmpHist)
            Kr(kk) = sum(tmpHist(1:kk))/sum(tmpHist);
        end
        Kr = Kr./lambda;
        Lr = sqrt(Kr./pi());
        Hr = Lr-Xr(1:end-1);
        RDFtable(ii,:) = [Hr(end), pXY(1), pXY(2), pXY(3)];
    else
        %no molicules in neighbourhood = H(r) undefined
        RDFtable(ii,:) = [NaN, pXY(1), pXY(2), pXY(3)];
    end
end

%% Generate RDFmap
%scale RDFtable; normalize to 1
minH = min(RDFtable(:,1));
tmpTable = RDFtable(:,1);
tmpTable = (tmpTable - minH) + 1;
tmpTable(isnan(tmpTable)) = 0;

%generate image
%scale
maxX = ceil(max(RDFtable(:,2)));
maxY = ceil(max(RDFtable(:,3)));
maxX = ceil(maxX/20); %scale to 20nm pixels
maxY = ceil(maxY/20); %scale to 20nm pixels

RDFmap = zeros(maxX, maxY); %empty grayscale image
for jj = 1:length(tmpTable)
    tPos(1) = ceil(RDFtable(jj,2)/20); %scale to 20nm pixels
    tPos(2) = ceil(RDFtable(jj,3)/20); %scale to 20nm pixels
    RDFmap(tPos(1),tPos(2)) = RDFmap(tPos(1),tPos(2),1) + tmpTable(jj);
end

maxH = max(RDFmap(:));
RDFmap = uint8((RDFmap./maxH).*255);

end %end main function