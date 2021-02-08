function [cropPos, cropTIFF, allPos, allTIFF] = LoadCrop (inputIm, cropRange, Conditions)

%==========================================================================
%                              FUNCTION
% Function to load between 1 and 3 position files, overlay the images and 
% produce a cropped and uncropped output files of both the positions & TIFs
%
%==========================================================================
%                             INPUT/OUTPUT
% Inputs:
%   -inputIm: cell array of strings of input file names; can have an
%             arbitrary number of cells, but only the first 3 entries will
%             be used in lists > 3 cells.  This field is mandatory.  Files
%             can be structured as:
%               -n rows of 2 columns (x/y positions)
%               -n rows of 3 columns (x/y/z positions)
%               -n rows of 4 columns (x/y/z positions & # photons)
%               -n rows of 5 columns (output of LeicaConv/qPALMconv)
%               -n rows of 5 columns (Leica GSD format)
%   -cropRange: area of image to keep, mandatory field, can be called as:
%       -0: no cropping (whole image)
%       -1: crop, using GUI
%       -[minX, maxX, minY, mxY]: min/max X/Y coordinates, in nm.  Must set
%                                 an appropriate scale in conditions, below
%   -Conditions: structure with the following fields.  If not called,
%                default settings are used.
%       -.ImScale: scale of TIFF images (in nm/pixel), default  = 20
%       -.saveTIFF: flag to save TIFF file (0 = no save),  default is 1.
%       -.saveUncropped: flag to save uncropped data, default is 0.
%       -.minPhotons: brightness filter; if non-zero, positions with fewer
%                     than minPhotons collected photons will be eliminated 
%                     from the data set.  Default  = 0;
%       -.minPrecission: same as min photons, but filters based on minimum
%                        precision rather than minimum number of collected
%                        photons. Default  = 0;
%       -.savePosFiles: Flag to automatically save position files. 0 = not
%                       saved. Default = 1;
%       -.saveTIFFs: Flag to automatically save TIFF's.  0 = not saved.
%                    Default  = 1
%       -.saveNonCrop: Flag to save non-cropped positions/tiffs;
%                      subsidiary to .savePosFiles and.saveTIFFs. Default =
%                      0
%       -.Ch1/2/3Label: Label for images in Ch1/2/3.  If left empty will
%                       be set as Ch1/2/3
%       -.fName: Output filename; will overwrite automated file name
%
% Notes: 
%   -File names are saved as Ch1_Ch2_Ch3.mat/.tif
%   -Input files must be pre-converted to .mat files
%
% Outputs:
%   -cropPos: structure with:
%       -.nCh: Number of channels in dataset (range: 1-3)
%       -.Ch(1/2/3): position matrices of each channel, each matrix is n x  
%                  5 matrix X/Y/Z coordinates,# of photons and precision. Matrix
%                  is of only those positions within the cropped region and 
%                  above the cutoff set by minPhotons/minPrecission, and is
%                  in nm units.
%       -.Conditions: same in input "Conditions"
%   -cropTIFF: 16-bit, 3 channel TIFF cropped image, built from .Ch1/2/3
%   -allPos: same as cropPos, but contains all positions that pass the
%            minPhoton/minPrecission filter.
%   -allTIFF: same as cropTIFF, but is of uncropped image
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

%% Check inputs & load files
if nargin ~= 2 && nargin ~=3
    edit 'GSDLoadCrop.m'
    error ('Input must be 2 or 3 variables.');
end

if ~iscell(inputIm)
    error ('First input must be a cell array of file names');
end

if length(cropRange) == 4
    cropRange(5:6) = 0;
end

if nargin == 2
    %Pan-function variables (keep the same across functions)
    Conditions.ImScale = 20; %nm/pixel in TIFF images
    Conditions.Ch1Label = 'Channel 1'; %channel name
    Conditions.Ch2Label = 'Channel 2'; %channel name
    Conditions.Ch3Label = 'Channel 3'; %channel name

    %LoadCrop-specific variables:
    Conditions.saveTIFF = 1; %save TIFF files
    Conditions.saveUncropped = 1; %save uncropped image
    Conditions.minPhotons = 0; %filter for minimal number of photons, 0 = no filter
    Conditions.minPrecission = 0; %filter for minimal precision (in nm), 0 = no filter
    Conditions.savePosFiles = 1; %save position (.mat) files
    Conditions.saveSingleCh = 0; %save single channels
end

nFiles = max(size(inputIm)); %determine number of images
if nFiles <= 3
    cropPos.nCh = nFiles;
    allPos.nCh = nFiles;
else
    cropPos.nCh = 3;
    allPos.nCh = 3;
    nFiles = 3;
end

cropPos.Conditions = Conditions;
allPos.Conditions = Conditions;

% load files & scale.
for ii=1:nFiles
    if exist(char(inputIm(ii))) == 2
        tPos = load(char(inputIm(ii)));
        tPos = tPos.outMat;
        if size(tPos,2) == 2
            allPos.Ch(ii).Pos(:,1:2) = tPos(:,1:2); %move x/y coords
            allPos.Ch(ii).Pos(:,3) = 0; %set z coordinates to zero
            allPos.Ch(ii).Pos(:,4) = 0; %set # photons to zero
            allPos.Ch(ii).Pos(:,5) = 0; %set precision to zero
        elseif size(tPos,2) == 3
            allPos.Ch(ii).Pos(:,1:3) = tPos(:,1:3); %move x/y/z coords
            allPos.Ch(ii).Pos(:,4) = 0; %set # photons to zero
            allPos.Ch(ii).Pos(:,5) = 0; %set precision to zero
        elseif size(tPos,2) == 4
            allPos.Ch(ii).Pos(:,:) = tPos(:,:); %move x/y/z coords & # photons
            allPos.Ch(ii).Pos(:,5) = 200./sqrt(allPos.Ch(ii).Pos(:,3)); %calculate precision
        elseif size(tPos,2) == 5
            allPos.Ch(ii).Pos = tPos; % copy over all data
        elseif size(tPos,2) == 9
            allPos.Ch(ii).Pos(:,1:2) = tPos(:,4:5); %x/y coordinates
            allPos.Ch(ii).Pos(:,3) = tPos(:,7); %z coordinates
            allPos.Ch(ii).Pos(:,4) = tPos(:,6); %copy-over number of photons
            allPos.Ch(ii).Pos(:,5) = 200./sqrt(allPos.Ch(ii).Pos(:,3)); %calculate precision
        else
            error ('Input matrix is not of a supported size.');
        end
    else
        error (['File ' char(inputIm(ii)) ' does not exist.']);
    end
end

clear tPos

%% Filter Inputs

%intensity filter data, if required
if Conditions.minPhotons ||Conditions.minPrecission
    if Conditions.minPrecission
        allPos.Ch(1).Pos(allPos.Ch(1).Pos(:,4)<Conditions.minPrecission,:) = []; %remove all entries below minPrecission
        allPos.Ch(2).Pos(allPos.Ch(2).Pos(:,4)<Conditions.minPrecission,:) = [];
        allPos.Ch(3).Pos(allPos.Ch(3).Pos(:,4)<Conditions.minPrecission,:) = [];
    else
        allPos.Ch(1).Pos(allPos.Ch(1).Pos(:,3)<Conditions.minPhotons,:) = []; %remove all entries below minPrecission
        allPos.Ch(2).Pos(allPos.Ch(2).Pos(:,3)<Conditions.minPhotons,:) = [];
        allPos.Ch(3).Pos(allPos.Ch(3).Pos(:,3)<Conditions.minPhotons,:) = [];
    end
end

%generate uncropped RGB image
for ii=1:nFiles
    maxX(ii) = max(allPos.Ch(ii).Pos(:,1));
    maxY(ii) = max(allPos.Ch(ii).Pos(:,2));
end
maxX = ceil((max(maxX))/Conditions.ImScale)+1; %set x-image size
maxY = ceil((max(maxY))/Conditions.ImScale)+1; %set y-image size

if nFiles == 1
    allTIFF = zeros(maxX, maxY, 1);
else
    allTIFF = zeros(maxX, maxY, 3);
end

for ii=1:nFiles
    tmpPos = ceil(allPos.Ch(ii).Pos(:,1:2)./Conditions.ImScale);
    tmpPos(tmpPos==0) = 1;
    for jj=1:length(tmpPos) 
        allTIFF(tmpPos(jj,1), tmpPos(jj,2),ii) = allTIFF(tmpPos(jj,1), tmpPos(jj,2),ii) + 1;
    end
    maxI(ii) = max(max(allTIFF(:,:,ii)));
    tmpPos = allTIFF(:,:,ii);
    tmpPos = tmpPos(:);
    tmpPos(tmpPos==0) = [];
    meanI(ii) = mean(tmpPos);
end

meanI = max(meanI);
if max(maxI)>65535 %if more than 16-bit
    for ii = 1:nFiles
        if maxI(ii)>65535
            allTIFF(:,:,ii) = (allTIFF(:,:,ii)./maxI(ii)).*65535; %rescale image intensity
        end
    end
end

%% Crop Image & Generate TIFF
if size (cropRange,2) == 6
    cropPos.cropRange = cropRange;
    for ii=1:nFiles
        cropPos.Ch(ii).Pos = allPos.Ch(ii).Pos(allPos.Ch(ii).Pos(:,1) >= cropRange(1) & allPos.Ch(ii).Pos(:,1) <= cropRange(2) & allPos.Ch(ii).Pos(:,2) >= cropRange(3) & allPos.Ch(ii).Pos(:,2) <= cropRange(4) & allPos.Ch(ii).Pos(:,3) >= cropRange(5) & allPos.Ch(ii).Pos(:,3) <= cropRange(6),:);
    end
elseif cropRange == 1 % use GUI to generate crop
    figure; imshow (allTIFF, [0 (meanI*8)]);
    Mask = imrect(gca);
    Mask = wait(Mask);
    close (gcf);
    Mask = double(Mask).*Conditions.ImScale;
    cropRange(3) = Mask(1);
    cropRange(4) = Mask(1) + Mask(3);
    cropRange(1) = Mask(2);
    cropRange(2) = Mask(2) + Mask(4);
    cropPos.cropRange = cropRange;
    for ii=1:nFiles
        cropPos.Ch(ii).Pos = allPos.Ch(ii).Pos(allPos.Ch(ii).Pos(:,1) >= cropRange(1) & allPos.Ch(ii).Pos(:,1) <= cropRange(2) & allPos.Ch(ii).Pos(:,2) >= cropRange(3) & allPos.Ch(ii).Pos(:,2) <= cropRange(4),:);
    end
else % no cropping
    cropPos.cropRange = cropRange;
    for ii=1:nFiles
        cropPos.Ch(ii).Pos = allPos.Ch(ii).Pos;
    end
end

%re-zero samples to origin + 1
for ii=1:nFiles
    minX(ii) = min(min(cropPos.Ch(ii).Pos(:,1)));
    minY(ii) = min(min(cropPos.Ch(ii).Pos(:,2)));
end
minX = min(minX)-1;
minY = min(minY)-1;

for ii=1:nFiles
    cropPos.Ch(ii).Pos(:,1) = cropPos.Ch(ii).Pos(:,1) - minX;
    cropPos.Ch(ii).Pos(:,2) = cropPos.Ch(ii).Pos(:,2) - minY;
end

%generate cropped TIFF image
for ii=1:nFiles
    maxX(ii) = max(cropPos.Ch(ii).Pos(:,1));
    maxY(ii) = max(cropPos.Ch(ii).Pos(:,2));
end
maxX = ceil((max(maxX))/Conditions.ImScale)+1; %set x-image size
maxY = ceil((max(maxY))/Conditions.ImScale)+1; %set y-image size

if nFiles == 1
    cropTIFF = zeros(maxX, maxY, 1);
else
    cropTIFF = zeros(maxX, maxY, 3);
end

for ii=1:nFiles
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

%% Save files & exit

allTIFF = uint16(allTIFF);

%generate output file name
if isfield (Conditions, 'fName')
    outName = Conditions.fName;
else
    if nFiles == 1
        outName = Conditions.Ch1Label;
    elseif nFiles == 2
        outName = [Conditions.Ch1Label '-' Conditions.Ch2Label];
    else
        outName = [Conditions.Ch1Label '-' Conditions.Ch2Label '-' Conditions.Ch3Label];
    end
end
    
if Conditions.saveTIFF
    imwrite (allTIFF, [outName '.tif'], 'TIFF');
    imwrite (cropTIFF, [outName ' - cropped.tif'], 'TIFF');
end

if Conditions.savePosFiles
    save ([outName ' - cropped.mat'], 'cropPos');
end

if Conditions.savePosFiles && Conditions.saveUncropped
    save ([outName ' - uncropped.mat'], 'allPos');
end

outName = {Conditions.Ch1Label; Conditions.Ch2Label; Conditions.Ch3Label};

if Conditions.saveSingleCh
    for ii=1:nFiles %save individual files
        outMat = cropPos.Ch(ii).Pos;
        save (['Single_Channel - ' char(outName(ii)) '.mat'], 'outMat');
    end
end

end %end function
