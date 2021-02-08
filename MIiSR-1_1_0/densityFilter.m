function densityFilter (CroppedPos, Conditions)

%==========================================================================
%                              FUNCTION
% This function to "de-speckle" GSD images by removing points which have 
% fewer than x * StDev neighbours relative to the image mean, looking in a 
% neighbourhood with a radius of Conditions.distFilt. This filtering is 
% useful for aiding DBSCAN and OPTICS analysis in the identification of 
% clusters within samples containing a large portion of unclustered 
% molecules.
%
%==========================================================================
%                              CAUTION
%   This script removes molecules located in regions of the image which
%   have low molecular density relative to the image average. Because this
%   filtering is not applied equally across the image, the use of this 
%   function should be strictly limited to pre-filtering of images prior to
%   DBSCAN and OPTICS segmentation. Data sets filtered with this function
%   should not be used for SAA, Cross-Correlation, Cross-Ripley's, or any 
%   other form quantitative analysis.
%
%==========================================================================
%                             INPUT/OUTPUT
% Inputs:
%   -CroppedPos = Position/intensity file produced by LoadCrop
%   -Conditions = Structure to set the varying analysis criteria:
%       -.distFilt: radius of local region or core distance used for 
%                   density calculations (in nm)
%       -.fMode: filter mode.  1 = filter on local density,  with the 
%                radius of the local area set t fistFilt (in nm) 2 = filter
%                on core distance, defined as the distance from each 
%                particle to the distFilt nearest molecule
%       -Ch1Filt = standard deviations below median density at which 
%               Channel 1 points are removed. [] = no filter
%       -.Ch2Filt = standard deviations below median density at which 
%               Channel 2 points are removed. [] = no filter
%       .Ch3Filt = standard deviations below median density at which 
%               Channel 3 points are removed. [] = no filter
%       -.saveTIFF = save TIFF-formatted filtered image. 0 = not saved, 1 =
%                    saved.
%
%  *** If Conditions is not called in the input, the following input values
%  are used (which will result in no filtering):

if (nargin == 1)
    display ('Using in-file conditions');
    Conditions.distFilt = 200; 
    Conditions.fMode = 1; 
    Conditions.Ch1Filt = [];
    Conditions.Ch2Filt = [];
    Conditions.Ch3Filt = [];
    Conditions.saveTIFF = 0;
end


% Outputs:
%   -Saves filtered .mat file, as ['Filtered - ' OrignoalFileName.mat] &
%    TIFF
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

%% Check Inputs
if nargin ~= 1 && nargin ~= 2
    open('GSDdensity.m');
    error ('Function must have 1 or 2 inputs');
end


if Conditions.fMode ~= 1 && Conditions.fMode ~= 2
    error ('Filter Mode must be either 1 (density) or 2 (core-distance)');
end

if Conditions.fMode == 2 && Conditions.distFilt < 1
    error ('Filtering in core-distance mode requires an integer entry > 1.');
end

if Conditions.fMode == 1
    Conditions.distFilt = ceil (Conditions.distFilt);
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

%load
load (CroppedPos);
if exist ('allPos', 'var')
    cropPos = allPos;
    clear allPos
    cropRange = [];
else
    cropRange = cropPos.cropRange;
end

nCh = cropPos.nCh;
tCond = cropPos.Conditions;

%load all channels (upto 3)
Ch1 = cropPos.Ch(1).Pos;
if nCh > 1
    Ch2 = cropPos.Ch(2).Pos;
end
if nCh == 3
    Ch3 = cropPos.Ch(3).Pos;
end

clear cropPos

%% Process for density
if ~isempty(Conditions.Ch1Filt)
    [cropPos.Ch(1).Pos, cropPos.MedDen(1)] =  filterForDenisty (Ch1, Conditions.fMode, Conditions.Ch1Filt, Conditions.distFilt);
else
    cropPos.Ch(1).Pos = Ch1;
end
clear Ch1

if nCh > 1 
    if ~isempty(Conditions.Ch2Filt)
        [cropPos.Ch(2).Pos, cropPos.MedDen(2)] =  filterForDenisty (Ch2, Conditions.fMode, Conditions.Ch2Filt, Conditions.distFilt);
    else
        cropPos.Ch(2).Pos = Ch2;
    end
    clear Ch2
end

if nCh == 3 
    if ~isempty(Conditions.Ch3Filt)
        [cropPos.Ch(3).Pos, cropPos.MedDen(3)] =  filterForDenisty (Ch3, Conditions.fMode, Conditions.Ch3Filt, Conditions.distFilt);
    else
        cropPos.Ch(3).Pos = Ch3;
    end
    clear Ch3
end

%% Save data to new output file
cropPos.nCh = nCh;
cropPos.cropRange = cropRange;
cropPos.Conditions = tCond;
cropPos.Conditions.distFilt = Conditions.distFilt;
cropPos.Conditions.fMode = Conditions.fMode; 
cropPos.Conditions.Ch1Filt = Conditions.Ch1Filt;
cropPos.Conditions.Ch2Filt = Conditions.Ch2Filt;
cropPos.Conditions.Ch3Filt = Conditions.Ch3Filt;
cropPos.Conditions.saveTIFF = Conditions.saveTIFF;

s = max(strfind(CroppedPos, '.'));
CroppedPos(s:end) = []; %remove extension

cropPos.OrigonalFile = CroppedPos;
save (['Filtered - ' CroppedPos], 'cropPos');

%save filtered tiff
if Conditions.saveTIFF 
    for ii=1:nCh
        maxX(ii) = max(cropPos.Ch(ii).Pos(:,1));
        maxY(ii) = max(cropPos.Ch(ii).Pos(:,2));
    end
    maxX = ceil((max(maxX))/Conditions.ImScale)+1; %set x-image size
    maxY = ceil((max(maxY))/Conditions.ImScale)+1; %set y-image size

    if nCh == 1
        cropTIFF = zeros(maxX, maxY, 1);
    else
        cropTIFF = zeros(maxX, maxY, 3);
    end

    for ii=1:nCh
        tmpPos = ceil(cropPos.Ch(ii).Pos(:,1:2)./Conditions.ImScale);
        for jj=1:length(tmpPos) 
            cropTIFF(tmpPos(jj,1), tmpPos(jj,2),ii) = cropTIFF(tmpPos(jj,1), tmpPos(jj,2),ii) + 1;
        end
        maxI(ii) = max(max(cropTIFF(:,:,ii)));
    end

    if max(maxI)>65535 %if more than 16-bit
        for ii = 1:nFiles
            if maxI(ii)>65535
                cropTIFF(:,:,ii) = (cropTIFF(:,:,ii)./maxI(ii)).*65535; %rescale image
            end
        end
    end

    cropTIFF = uint16(cropTIFF);

    imwrite (cropTIFF, ['Filtered - ' CroppedPos], 'TIFF');
end %end save TIFF

%% Clean-up
if openPool == 1 %close pool of opened by function
    delete (gcp);
end

end %end main function


%% Sub-Functions
function [outMat, MedDen] = filterForDenisty (inMat, fMode, StDev, distFilt)
    m = size(inMat,1);
    nPart = zeros(m,1); %number of particles

    if fMode == 1
        parfor ii=1:m
            tPos = inMat(:,1:3);
            tPos(ii,:) = []; %eliminate self-position
            D = sqrt((tPos(:,1)-inMat(ii,1)).^2 + (tPos(:,2)-inMat(ii,2)).^2 + (tPos(:,3)-inMat(ii,3)).^2);
            nPart(ii) = sum(D<=distFilt);
        end

    else
        parfor ii=1:m
            tPos = inMat(:,1:3);
            tPos(ii,:) = []; %eliminate self-position
            D = sqrt((tPos(:,1)-inMat(ii,1)).^2 + (tPos(:,2)-inMat(ii,2)).^2 + (tPos(:,3)-inMat(ii,3)).^2);
            nPart(ii) = D(distFilt);
        end
    end

    % Filter Data
    MedDen = median(nPart);
    dSTD = std(nPart);
    cutOff = MedDen + (dSTD*StDev);

    if fMode == 1
        outMat = inMat(nPart >= cutOff,:);
    else
        outMat = inMat(nPart <= cutOff,:);
    end

    if fMode == 1
        MedDen = MedDen/(pi()*distFilt^2);
    else
        MedDen = NaN;
    end

end %end function
