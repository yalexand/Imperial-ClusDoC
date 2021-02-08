function SAAstruct = SAA3col(CroppedPos, nCores, Conditions)
%==========================================================================
%                              FUNCTION
% Function to perform spatial association analysis on GSD position files
% using a position/intensity file produced from GSDLoadCrop. 
%
% Function which performs a spatial association analysis of
% single-molecules images.  The Function measures the distances between the 
% particles in Im1 and the closest particle in Im2 & Im3, and then compares 
% the histogram of measured distances to a histogram derived from 
% randomized particle positions within an image of the same area as the 
% input.
%
%==========================================================================
%                             INPUT/OUTPUT
% Inputs:
%   -CroppedPos = Position/intensity file produced by GSDLoadCrop
%   -nCores = number of cores to use for processing data.  Set at 0 to use
%   default settings, set at 1 to prevent parallel processing.
%   -Conditions = Structure to set the varying analysis criteria:
%       -Iterations = number of randomized images to generate to form
%        randomized distance histogram
%       -Scale = image scale (nm/pixel, um/pixel, etc).
%       -SD = standard deviations for SMSA algorithm.  p = 1.65 for 90% 
%             cut-off, p = 2.00 for 95% cut-off
%       -Cutoff = value to add to the cutoff determined by function.  Used
%        to account for known fluorophore spacing, deal with image.  Must
%        be in same units as Conditions.Scale.
%        registration defects, etc.  Set to 0 if not needed.
%
%  *** If Conditions is not called in the input, the following input values
%  are used:

if (nargin == 2)
    display ('Using in-file conditions');

    %Pan-function variables (keep the same across functions)
    Conditions.Ch1Label = 'Channel 1'; %label in image 1
    Conditions.Ch2Label = 'Channel 2'; %label in image 2
    Conditions.Ch3Label = 'Channel 3'; %label in image 3

    %SAA-specific variables
    Conditions.AnalyzedFraction = 1; %Fraction of points to use in quantification% set at less than 1 for more rapid, but less accurate quantification
    Conditions.SD = 1.65; %Standard deviation of error.  1.65 = p of 0.90; 2.00 = p of 0.95.
    Conditions.Cutoff = 0; %Value to add to SAA, for known separations or poor image registration
    Conditions.BinMax = 1000; %maximum distance to assess
    Conditions.Iterations = 1; %randomization iterations for Simulated Random Positions calculations
    Conditions.gFilter = 3; %Diameter of smoothing filter for contour graphs.
    Conditions.ContourDensity = 25; %Contour line density on plot
    Conditions.polyRand = 0; %Perform randomization over an area determined using convhull; used to ensure that clear 
                             % areas created by invisible objects are not
                             % included in SAA calculations. 0 = off; 1 = on
end

% Outputs:
%
%   -OutputData
%       .Im1MeaPhotons/Im2/Im3: Average photons in Im1/2/3
%       .Im1MeanPrecission/Im2/Im3: Mean precision of Im1/2/3.
%           -Calculated as 200/sqrt(number photons)
%       .sigmaRMS: Root-mean squared of Im1/Im2 precision
%       .SAAcutoff: cut-off for SAA analysis, 
%       .distTabelm1/Im2: nearest neighbour distances for Im1/Im2
%       .Im1Graph/Im2: Histogram of nearest neighbour distances for Im1/Im2
%       .randTablIm1/Im2: nearest neighbour table for randomized data.  One
%        column per Conditions.Repeats (iterations).  Last column is
%        average of iterations.
%       .randTableIm1/Im2: Histogram of nearest neighbour distances for
%        randomized data sets
%       .xAxisLabels: x-axis values for histograms
%       .BarData: particle counts for bar graph (figure 3)
%       .eBarData: data for enrichment bar graph (figure 4)
%       .fBarData: fraction of particles that are co-localized (figure 5)
%
%   -RawData
%       .Im1/Im2/Im3 - Image 1/2/3 file name
%       .Conditions - Analysis conditions used
%       .Im1Pos/Im2Pos - *Post-filtering* particle positions used for
%        analysis
%
%   Note: RawData & OutputData are automatically saved as a .m file.
%
%   Figures:
%   -Figure 1: SAA plot for Im1:Im2
%   -Figure 2: SAA plot for Im2:Im1
%   -Figure 3: Bar graph of number/Fraction of particles below SAA
%   -Figure 4: Fold-enrichment of particles below SAA relative to
%              randomized populations
%   -Figure 5: Figures 1-4 compiled into a 2x2 figure, for printing.
%              Outputted as both a .fig and .pdf
%
%==========================================================================
%                       NOTES & DEPENDENCIES
%   Dependencies: A recent copy of Matlab (2010a or newer) plus the 
%                 parallel processing toolbox. This code is processor-
%                 intensive, thus a high-end, multi-processor workstation 
%                 is recommended.
%
%  1) The cutoff used to measure co-localized fraction is equal too the
%     sum of the root-mean-squared of the measured precision for each
%     channel, plus Conditions.Cutoff.
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
ttime = tic;

if (nargin ~= 2) && (nargin ~= 3)
    open ('GSD_SAA3d.m');
    error ('Input is either 2 or 3 variables.  See GSD_SAA2d.m for help');
end

%check for files

if ~(exist (CroppedPos, 'file'))
    error ('Input file is not in working directory');
end

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

%% Load Data & Prepare for Processing
RawData.Conditions = Conditions;
xAxisLabels = [0:Conditions.BinMax/100:Conditions.BinMax];

load (CroppedPos);
if exist ('cropPos', 'var')
    RawData.Im1Pos = cropPos.Ch(1).Pos;
    RawData.Im2Pos = cropPos.Ch(2).Pos;
    RawData.Im3Pos = cropPos.Ch(3).Pos;
else
    RawData.Im1Pos = allPos.Ch(1).Pos;
    RawData.Im2Pos = allPos.Ch(2).Pos;
    RawData.Im3Pos = allPos.Ch(3).Pos;
end
clear cropPos


%reduce data to subset, if required
if Conditions.AnalyzedFraction < 1
    aT = tic;
    display (' ');
    display ('Reducing data size.');
    PosKeep = rand(length(RawData.Im1Pos),1);
    RawData.Im1Pos = RawData.Im1Pos(PosKeep <= Conditions.AnalyzedFraction,:);
    PosKeep = rand(length(RawData.Im2Pos),1);
    RawData.Im2Pos = RawData.Im2Pos(PosKeep <= Conditions.AnalyzedFraction,:);
    PosKeep = rand(length(RawData.Im3Pos),1);
    RawData.Im3Pos = RawData.Im3Pos(PosKeep <= Conditions.AnalyzedFraction,:);
    display (['. . .Data reduction completed in ' num2str(toc(aT)) ' seconds.']);
    clear PosKeep aT;
end %end subset generation


%determine if data is 2D or 3D
nDims = max([max(RawData.Im1Pos(:,3)), max(RawData.Im2Pos(:,3)), max(RawData.Im3Pos(:,3))]);
if nDims ~= 0
    nDims = 3;
else
    nDims = 2;
end

%Calculate precission & sigma-RMS
OutputData.Im1MeanPhotons = mean(RawData.Im1Pos(:,4));
OutputData.Im1MeanPrecission = mean(RawData.Im1Pos(:,5));
OutputData.Im2MeanPhotons = mean(RawData.Im2Pos(:,4));
OutputData.Im2MeanPrecission = mean(RawData.Im2Pos(:,5));
OutputData.Im3MeanPhotons = mean(RawData.Im3Pos(:,4));
OutputData.Im3MeanPrecission = mean(RawData.Im3Pos(:,5));
OutputData.sigmaRMS3Ch = sqrt((OutputData.Im1MeanPrecission^2) + (OutputData.Im2MeanPrecission^2) + (OutputData.Im3MeanPrecission^2));
OutputData.sigmaRMSCh1_2 = sqrt((OutputData.Im1MeanPrecission^2) + (OutputData.Im2MeanPrecission^2));
OutputData.sigmaRMSCh1_3 = sqrt((OutputData.Im1MeanPrecission^2) + (OutputData.Im3MeanPrecission^2));
OutputData.sigmaRMSCh2_3 = sqrt((OutputData.Im2MeanPrecission^2) + (OutputData.Im3MeanPrecission^2));
OutputData.SAAcutoff3Ch = (OutputData.sigmaRMS3Ch*Conditions.SD) + Conditions.Cutoff;
OutputData.SAAcutoffCh1_2 = (OutputData.sigmaRMSCh1_2*Conditions.SD) + Conditions.Cutoff;
OutputData.SAAcutoffCh1_3 = (OutputData.sigmaRMSCh1_3*Conditions.SD) + Conditions.Cutoff;
OutputData.SAAcutoffCh2_3 = (OutputData.sigmaRMSCh2_3*Conditions.SD) + Conditions.Cutoff;

%% Generate Distance Tables
% Using fastest eucld. distance equation - http://www.mathworks.com/matlabcentral/fileexchange/71-distance-m
display ('...Performing 2-Channel SAA Analysis ');
display ('   ...Generating Distance Tables.');

display (' ');
display ('Generating Distance Tables.');
nPartIm1 = size(RawData.Im1Pos,1);
nPartIm2 = size(RawData.Im2Pos,1);
nPartIm3 = size(RawData.Im3Pos,1);
Range = Conditions.BinMax;

se1 = strel ('disk', Conditions.gFilter);
se1 = se1.getnhood;
se1 = se1./sum(se1(:)); %structuring element for smoothing of contour plot

%pre-allocate for speed
distTableIm1 = zeros(nPartIm1,2);
distTableIm2 = zeros(nPartIm2,2);
distTableIm3 = zeros(nPartIm3,2);
Pos1 = RawData.Im1Pos(:,1:3)';
Pos2 = RawData.Im2Pos(:,1:3)';
Pos3 = RawData.Im3Pos(:,1:3)';

%================== Process Pos1 ================
aT = tic;
parfor f=1:nPartIm1
    [a, b] =  nearestNeighbour(Pos1(:,f), Pos2, Pos3, Range);
    distTableIm1(f,:) = [a, b];
end
%Data for Graphs
OutputData.distTableIm1 = distTableIm1';
%2-channel histograms
GraphData.Im1Graph1_2 = hist(distTableIm1(:,1)', xAxisLabels);
GraphData.Im1Graph1_2 = GraphData.Im1Graph1_2./nPartIm1;
GraphData.Im1Graph1_3 = hist(distTableIm1(:,2)', xAxisLabels);
GraphData.Im1Graph1_3 = GraphData.Im1Graph1_3./nPartIm1;

%contour plots
tPlot = zeros(Range, Range);
for i=1:length (distTableIm1)
    pX(1) = uint16(distTableIm1(i,1))+1;
    pX(2) = uint16(distTableIm1(i,2))+1;
    if (pX(1) <= Range) && (pX(2) <= Range)
        tPlot(pX(1),pX(2)) = tPlot(pX(1),pX(2)) + 1;
    end
end
GraphData.Im1Contour = tPlot./Conditions.Iterations;
GraphData.Im1ContourFiltered = imfilter(tPlot,se1); 
GraphData.Im1Contour = GraphData.Im1Contour ./nPartIm1;
GraphData.Im1ContourFiltered  = GraphData.Im1ContourFiltered ./nPartIm1;


display ('      -Distance table calculated for dataset 1.');
display (['          -completed in ' num2str(toc(aT)) ' seconds.']);
clear aT a b tPlot distTableIm1 pX tPlot


%================== Process Pos2 ================
aT = tic;
parfor f=1:nPartIm2
    [a, b] =  nearestNeighbour(Pos2(:,f), Pos1, Pos3, Range);
    distTableIm2(f,:) = [a, b];
end
%Data for Graphs
OutputData.distTableIm2 = distTableIm2';
%2-channel histograms
GraphData.Im2Graph2_1 = hist(distTableIm2(:,1)', xAxisLabels);
GraphData.Im2Graph2_1 = GraphData.Im2Graph2_1./nPartIm2;
GraphData.Im2Graph2_3 = hist(distTableIm2(:,2)', xAxisLabels);
GraphData.Im2Graph2_3 = GraphData.Im2Graph2_3./nPartIm2;

%contour plots
tPlot = zeros(Range, Range);
for i=1:length (distTableIm2)
    pX(1) = uint16(distTableIm2(i,1))+1;
    pX(2) = uint16(distTableIm2(i,2))+1;
    if (pX(1) <= Range) && (pX(2) <= Range)
        tPlot(pX(1),pX(2)) = tPlot(pX(1),pX(2)) + 1;
    end
end
GraphData.Im2Contour = tPlot./Conditions.Iterations;
GraphData.Im2ContourFiltered = imfilter(tPlot,se1); 
GraphData.Im2Contour = GraphData.Im2Contour ./nPartIm2;
GraphData.Im2ContourFiltered  = GraphData.Im2ContourFiltered ./nPartIm2;

display ('      -Distance table calculated for dataset 2.');
display (['          -completed in ' num2str(toc(aT)) ' seconds.']);
clear aT a b tPlot distTableIm2 pX tPlot


%================== Process Pos3 ================
aT = tic;
parfor f=1:nPartIm3
    [a, b] =  nearestNeighbour(Pos3(:,f), Pos1, Pos2, Range);
    distTableIm3(f,:) = [a, b];
end
%Data for Graphs
OutputData.distTableIm3 = distTableIm3';
%2-OutputData.SAAcutoffCh1_2channel histograms
GraphData.Im3Graph3_1 = hist(distTableIm3(:,1)', xAxisLabels);
GraphData.Im3Graph3_1 = GraphData.Im3Graph3_1./nPartIm3;
GraphData.Im3Graph3_2 = hist(distTableIm3(:,2)', xAxisLabels);
GraphData.Im3Graph3_2 = GraphData.Im3Graph3_2./nPartIm3;

%contour plots
tPlot = zeros(Range, Range);
for i=1:length (distTableIm3)
    pX(1) = uint16(distTableIm3(i,1))+1;
    pX(2) = uint16(distTableIm3(i,2))+1;
    if (pX(1) <= Range) && (pX(2) <= Range)
        tPlot(pX(1),pX(2)) = tPlot(pX(1),pX(2)) + 1;
    end
end
GraphData.Im3Contour = tPlot./Conditions.Iterations;
GraphData.Im3ContourFiltered = imfilter(tPlot,se1); 
GraphData.Im3Contour = GraphData.Im3Contour ./nPartIm3;
GraphData.Im3ContourFiltered  = GraphData.Im3ContourFiltered ./nPartIm3;

display ('      -Distance table calculated for dataset 3.');
display (['          -completed in ' num2str(toc(aT)) ' seconds.']);
clear aT a b tPlot distTableIm3 pX tPlot


%% Random Populations
if Conditions.Iterations
    display (' ');
    display ('   ...Generating Randomized Distance Tables.');
    
    % ========= random tables =============
    % determine image area to use for randomization
    if Conditions.polyRand % if ROI is refined using convenx hull
        if nDims == 2
            tPos = vertcat(RawData.Im1Pos(:,1:2), RawData.Im2Pos(:,1:2), RawData.Im3Pos(:,1:2));
            [~, tPos] = convhull(tPos); %area of convex hull
            dX = sqrt(tPos);
            dY = dX;
            dZ = 0;
            clear tPos
        else % if using user-set rectangular area
            tPos = vertcat(RawData.Im1Pos(:,1:3), RawData.Im2Pos(:,1:3), RawData.Im3Pos(:,1:3));
            [~, tPos] = convhull(tPos); %volume of conex hull
            dX = nthroot(tPos,3);
            dY = dX;
            dZ = dX;
            clear tPos
        end
    else
        dX = max([max(Pos1(1,:)), max(Pos2(1,:))]) - min([min(Pos1(1,:)), min(Pos2(1,:))]);
        dY = max([max(Pos1(2,:)), max(Pos2(2,:))]) - min([min(Pos1(2,:)), min(Pos2(2,:))]);
        dZ = max([max(Pos1(3,:)), max(Pos2(3,:))]) - min([min(Pos1(3,:)), min(Pos2(3,:))]);
    end

    RandSize1 = size(Pos1);
    RandSize2 = size(Pos2);
	RandSize3 = size(Pos3);


    for i=1:Conditions.Iterations
        % prepare random image for channel 1
        randTable1(i).a = rand(RandSize1);
        randTable1(i).a(1,:) = randTable1(i).a(1,:).*dX;
        randTable1(i).a(2,:) = randTable1(i).a(2,:).*dY;
        randTable1(i).a(3,:) = randTable1(i).a(3,:).*dZ;
        
        % prepare random image for channel 2
        randTable2(i).a = rand(RandSize2);
        randTable2(i).a(1,:) = randTable2(i).a(1,:).*dX;
        randTable2(i).a(2,:) = randTable2(i).a(2,:).*dY;
        randTable2(i).a(3,:) = randTable2(i).a(3,:).*dZ;
        
        % prepare random image for channel 3
        randTable3(i).a = rand(RandSize3);
        randTable3(i).a(1,:) = randTable3(i).a(1,:).*dX;
        randTable3(i).a(2,:) = randTable3(i).a(2,:).*dY;
        randTable3(i).a(3,:) = randTable3(i).a(3,:).*dZ;
    end
    
     clear dX dY dZ RandSize*;
    
    % == Process Im1 ==
    aT = tic;
    rTable1 = zeros(nPartIm1,Conditions.Iterations);
    rTable2 = zeros(nPartIm1,Conditions.Iterations);
    for i=1:Conditions.Iterations
        display (['     -' Conditions.Ch1Label ', iteration ' num2str(i) ' of ' num2str(Conditions.Iterations) '.']);
        tTable1 = randTable1(i).a; 
        tTable2 = randTable2(i).a;
        tTable3 = randTable3(i).a;
        parfor f=1:nPartIm1
            [a, b] =  nearestNeighbour(Pos1(:,f), tTable2, tTable3, Range);
            rTable1(f,i) = a;
            rTable2(f,i) = b;
        end
    end
    %data for graphs
    Randomizations.randTableIm1_2 = rTable1';
    Randomizations.randTableIm1_3 = rTable2';
    Randomizations.RandGraph1_2 = hist(rTable1, xAxisLabels);
    Randomizations.RandGraph1_3 = hist(rTable2, xAxisLabels);
    Randomizations.RandGraph1_2 = Randomizations.RandGraph1_2(:,i)./nPartIm1;
    Randomizations.RandGraph1_3 = Randomizations.RandGraph1_3(:,i)./nPartIm1;
    Randomizations.RandGraph1_2(:,Conditions.Iterations+1) = mean(Randomizations.RandGraph1_2, 2);
    Randomizations.RandGraph1_3(:,Conditions.Iterations+1) = mean(Randomizations.RandGraph1_3, 2);
    
    %contour plot
    tPlot = zeros(Range, Range);
    for j=1:Conditions.Iterations
        for i=1:nPartIm1
            pX(1) = uint16(rTable1(i,j))+1;
            pX(2) = uint16(rTable2(i,j))+1;
            if (pX(1) <= Range) && (pX(2) <= Range)
                tPlot(pX(1),pX(2)) = tPlot(pX(1),pX(2)) + 1;
            end
        end
    end
    Randomizations.RandGraph1Contoured = tPlot./Conditions.Iterations;
    Randomizations.RandGraph1ContouredFiltered = imfilter(tPlot,se1); 
    Randomizations.RandGraph1Contoured = Randomizations.RandGraph1Contoured./nPartIm1;
    Randomizations.RandGraph1ContouredFiltered  = Randomizations.RandGraph1ContouredFiltered ./nPartIm1;

    clear randTableIm1* tTable* tPlot pX rTable*
    display (['          -' Conditions.Ch1Label ' randomization iterations completed in ' num2str(toc(aT)) ' seconds.']);
    clear aT
    
    % == Process Im2 ==
    aT = tic;
    rTable1 = zeros(nPartIm2,Conditions.Iterations);
    rTable2 = zeros(nPartIm2,Conditions.Iterations);
    for i=1:Conditions.Iterations
        display (['     -' Conditions.Ch2Label ', iteration ' num2str(i) ' of ' num2str(Conditions.Iterations) '.']);
        tTable1 = randTable1(i).a; 
        tTable2 = randTable2(i).a;
        tTable3 = randTable3(i).a;
        parfor f=1:nPartIm2
            [a, b] =  nearestNeighbour(Pos2(:,f), tTable1, tTable3, Range);
            rTable1(f,i) = a;
            rTable2(f,i) = b;
        end
    end
    %data for graphs
    Randomizations.randTableIm2_1 = rTable1';
    Randomizations.randTableIm2_3 = rTable2';
    Randomizations.RandGraph2_1 = hist(rTable1, xAxisLabels);
    Randomizations.RandGraph2_3 = hist(rTable2, xAxisLabels);
    Randomizations.RandGraph2_1 = Randomizations.RandGraph2_1(:,i)./nPartIm2;
    Randomizations.RandGraph2_3 = Randomizations.RandGraph2_3(:,i)./nPartIm2;
    Randomizations.RandGraph2_1(:,Conditions.Iterations+1) = mean(Randomizations.RandGraph2_1, 2);
    Randomizations.RandGraph2_3(:,Conditions.Iterations+1) = mean(Randomizations.RandGraph2_3, 2);
    
    %contour plot
    tPlot = zeros(Range, Range);
    for j=1:Conditions.Iterations
        for i=1:nPartIm2
            pX(1) = uint16(rTable1(i,j))+1;
            pX(2) = uint16(rTable2(i,j))+1;
            if (pX(1) <= Range) && (pX(2) <= Range)
                tPlot(pX(1),pX(2)) = tPlot(pX(1),pX(2)) + 1;
            end
        end
    end
    Randomizations.RandGraph2Contoured = tPlot./Conditions.Iterations;
    Randomizations.RandGraph2ContouredFiltered = imfilter(tPlot,se1); 
    Randomizations.RandGraph2Contoured = Randomizations.RandGraph2Contoured./nPartIm2;
    Randomizations.RandGraph2ContouredFiltered  = Randomizations.RandGraph2ContouredFiltered ./nPartIm2;

    clear randTableIm1* tTable* tPlot pX rTable*
    display (['          -' Conditions.Ch2Label ' randomization iterations completed in ' num2str(toc(aT)) ' seconds.']);
    clear aT
    
    % == Process Im3 ==
    aT = tic;
    rTable1 = zeros(nPartIm3,Conditions.Iterations);
    rTable2 = zeros(nPartIm3,Conditions.Iterations);
    for i=1:Conditions.Iterations
        display (['     -' Conditions.Ch3Label ', iteration ' num2str(i) ' of ' num2str(Conditions.Iterations) '.']);
        tTable1 = randTable1(i).a; 
        tTable2 = randTable2(i).a;
        tTable3 = randTable3(i).a;
        parfor f=1:nPartIm3
            [a, b] =  nearestNeighbour(Pos3(:,f), tTable1, tTable2, Range);
            rTable1(f,i) = a;
            rTable2(f,i) = b;
        end
    end
    %data for graphs
    Randomizations.randTableIm3_1 = rTable1';
    Randomizations.randTableIm3_2 = rTable2';
    Randomizations.RandGraph3_1 = hist(rTable1, xAxisLabels);
    Randomizations.RandGraph3_2 = hist(rTable2, xAxisLabels);
    Randomizations.RandGraph3_1 = Randomizations.RandGraph3_1(:,i)./nPartIm3;
    Randomizations.RandGraph3_2 = Randomizations.RandGraph3_2(:,i)./nPartIm3;
    Randomizations.RandGraph3_1(:,Conditions.Iterations+1) = mean(Randomizations.RandGraph3_1, 2);
    Randomizations.RandGraph3_2(:,Conditions.Iterations+1) = mean(Randomizations.RandGraph3_2, 2);
    
    %contour plot
    tPlot = zeros(Range, Range);
    for j=1:Conditions.Iterations
        for i=1:nPartIm3
            pX(1) = uint16(rTable1(i,j))+1;
            pX(2) = uint16(rTable2(i,j))+1;
            if (pX(1) <= Range) && (pX(2) <= Range)
                tPlot(pX(1),pX(2)) = tPlot(pX(1),pX(2)) + 1;
            end
        end
    end
    Randomizations.RandGraph3Contoured = tPlot./Conditions.Iterations;
    Randomizations.RandGraph3ContouredFiltered = imfilter(tPlot,se1); 
    Randomizations.RandGraph3Contoured = Randomizations.RandGraph3Contoured./nPartIm3;
    Randomizations.RandGraph3ContouredFiltered  = Randomizations.RandGraph3ContouredFiltered ./nPartIm3;

    clear randTableIm1* tTable* tPlot pX rTable*
    display (['          -' Conditions.Ch3Label ' randomization iterations completed in ' num2str(toc(aT)) ' seconds.']);
    clear aT
end %end if randomizations

clear Pos1 Pos2 Pos3 se1

%% Graphs
OutputData.xAxisLabels = xAxisLabels;

%Plot Im1:Im2/3
YMax(1,1:2) = [max(GraphData.Im1Graph1_2), max(GraphData.Im1Graph1_3)];
YMax(2,1:2) = [max(Randomizations.RandGraph1_2(:,end)), max(Randomizations.RandGraph1_3(:,end))];
YMax = max(YMax);
YMax = ceil(YMax.*20)/20; %round up to nearest 0.05

h = figure(1);
tPos = get(figure(1), 'Position');
set(figure(1), 'Position', [75, tPos(2), 1120, 1450]);
set (figure(1), 'Color', 'w');
text(0.35,1.05,['Two Channel Colocalizations: ' Conditions.Ch1Label ', ' Conditions.Ch2Label ' & ' Conditions.Ch3Label]);

subplot (3,2,1);
plot (OutputData.xAxisLabels, GraphData.Im1Graph1_2,'k');
hold on;
plot (OutputData.xAxisLabels, Randomizations.RandGraph1_2(:,end), '--k');
line ([OutputData.sigmaRMSCh1_2, OutputData.sigmaRMSCh1_2], [0, YMax(1)], 'LineStyle', '-.', 'Color', 'k');
xlabel ([Conditions.Ch1Label ':' Conditions.Ch2Label ' Separation Distance']);
ylabel ('Fraction of Particles');
axis ([0, Conditions.BinMax, 0, YMax(1)]);
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
hold off;

subplot (3,2,2);
plot (OutputData.xAxisLabels, GraphData.Im1Graph1_3,'k');
hold on;
plot (OutputData.xAxisLabels, Randomizations.RandGraph1_3(:,end), '--k');
line ([OutputData.sigmaRMSCh1_3, OutputData.sigmaRMSCh1_3], [0, YMax(2)], 'LineStyle', '-.', 'Color', 'k');
xlabel ([Conditions.Ch1Label ':' Conditions.Ch3Label ' Separation Distance']);
ylabel ('Fraction of Particles');
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
axis ([0, Conditions.BinMax, 0, YMax(2)]);
hold off;


%Plot Im2:Im1/3
YMax(1,1:2) = [max(GraphData.Im2Graph2_1), max(GraphData.Im2Graph2_3)];
YMax(2,1:2) = [max(Randomizations.RandGraph2_1(:,end)), max(Randomizations.RandGraph2_3(:,end))];
YMax = max(YMax);
YMax = ceil(YMax.*20)/20; %round up to nearest 0.05

subplot (3,2,3);
plot (OutputData.xAxisLabels, GraphData.Im2Graph2_1,'k');
hold on;
plot (OutputData.xAxisLabels, Randomizations.RandGraph2_1(:,end), '--k');
line ([OutputData.sigmaRMSCh1_2, OutputData.sigmaRMSCh1_2], [0, YMax(1)], 'LineStyle', '-.', 'Color', 'k');
xlabel ([Conditions.Ch2Label ':' Conditions.Ch1Label ' Separation Distance']);
ylabel ('Fraction of Particles');
axis ([0, Conditions.BinMax, 0, YMax(1)]);
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
hold off;

subplot (3,2,4);
plot (OutputData.xAxisLabels, GraphData.Im2Graph2_1,'k');
hold on;
plot (OutputData.xAxisLabels, Randomizations.RandGraph2_1(:,end), '--k');
line ([OutputData.sigmaRMSCh2_3, OutputData.sigmaRMSCh2_3], [0, YMax(2)], 'LineStyle', '-.', 'Color', 'k');
xlabel ([Conditions.Ch2Label ':' Conditions.Ch3Label ' Separation Distance']);
ylabel ('Fraction of Particles');
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
axis ([0, Conditions.BinMax, 0, YMax(2)]);
hold off;

%Plot Im3:Im1/2
YMax(1,1:2) = [max(GraphData.Im3Graph3_1), max(GraphData.Im3Graph3_2)];
YMax(2,1:2) = [max(Randomizations.RandGraph3_1(:,end)), max(Randomizations.RandGraph3_2(:,end))];
YMax = max(YMax);
YMax = ceil(YMax.*20)/20; %round up to nearest 0.05

subplot (3,2,5);
plot (OutputData.xAxisLabels, GraphData.Im3Graph3_1,'k');
hold on;
plot (OutputData.xAxisLabels, Randomizations.RandGraph3_1(:,end), '--k');
line ([OutputData.sigmaRMSCh1_3, OutputData.sigmaRMSCh1_3], [0, YMax(1)], 'LineStyle', '-.', 'Color', 'k');
xlabel ([Conditions.Ch3Label ':' Conditions.Ch1Label ' Separation Distance']);
ylabel ('Fraction of Particles');
axis ([0, Conditions.BinMax, 0, YMax(1)]);
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
hold off;

subplot (3,2,6);
plot (OutputData.xAxisLabels, GraphData.Im3Graph3_2,'k');
hold on;
plot (OutputData.xAxisLabels, Randomizations.RandGraph3_2(:,end), '--k');
line ([OutputData.sigmaRMSCh2_3, OutputData.sigmaRMSCh2_3], [0, YMax(2)], 'LineStyle', '-.', 'Color', 'k');
xlabel ([Conditions.Ch3Label ':' Conditions.Ch2Label ' Separation Distance']);
ylabel ('Fraction of Particles');
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
axis ([0, Conditions.BinMax, 0, YMax(2)]);
hold off;

text(-600,1.4,['Two Channel Colocalizations: ' Conditions.Ch1Label ', ' Conditions.Ch2Label ' & ' Conditions.Ch3Label]);
saveas (figure(1), 'SAA - Two-Channel Analysis');
set (figure(1), 'PaperType', 'usletter')
save2pdf ('SAA - Two-Channel Analysis.pdf', figure(1), 600);
close (h);

%bar data 
%  -arranged [triple, Im1:Im2, Im1:Im3] 
%  -2nd line is randomized data

rT1(1,:) = Randomizations.randTableIm1_2(:);
rT1(2,:) = Randomizations.randTableIm1_3(:);

GraphData.BarDataIm1(1,1) = sum((OutputData.distTableIm1(1,:) <= OutputData.SAAcutoff3Ch) & (OutputData.distTableIm1(2,:) <= OutputData.SAAcutoff3Ch))/nPartIm1;
GraphData.BarDataIm1(1,2) = sum((OutputData.distTableIm1(1,:) <= OutputData.SAAcutoff3Ch) & (OutputData.distTableIm1(2,:) > OutputData.SAAcutoff3Ch))/nPartIm1;
GraphData.BarDataIm1(1,3) = sum((OutputData.distTableIm1(1,:) > OutputData.SAAcutoff3Ch) & (OutputData.distTableIm1(2,:) <= OutputData.SAAcutoff3Ch))/nPartIm1;
GraphData.BarDataIm1(1,4) = 1-sum(GraphData.BarDataIm1(1,1:3));

GraphData.BarDataIm1(2,1) = sum((rT1(1,:) <= OutputData.SAAcutoff3Ch) & (rT1(2,:) <= OutputData.SAAcutoff3Ch))/(Conditions.Iterations*nPartIm1);
GraphData.BarDataIm1(2,2) = sum((rT1(1,:) <= OutputData.SAAcutoff3Ch) & (rT1(2,:) > OutputData.SAAcutoff3Ch))/(Conditions.Iterations*nPartIm1);
GraphData.BarDataIm1(2,3) = sum((rT1(1,:) > OutputData.SAAcutoff3Ch) & (rT1(2,:) <= OutputData.SAAcutoff3Ch))/(Conditions.Iterations*nPartIm1);
GraphData.BarDataIm1(2,4) = 1-sum(GraphData.BarDataIm1(2,1:3));
clear rT1

rT2(1,:) = Randomizations.randTableIm2_1(:);
rT2(2,:) = Randomizations.randTableIm2_3(:);

GraphData.BarDataIm2(1,1) = sum((OutputData.distTableIm2(1,:) <= OutputData.SAAcutoff3Ch) & (OutputData.distTableIm2(2,:) <= OutputData.SAAcutoff3Ch))/nPartIm2;
GraphData.BarDataIm2(1,2) = sum((OutputData.distTableIm2(1,:) <= OutputData.SAAcutoff3Ch) & (OutputData.distTableIm2(2,:) > OutputData.SAAcutoff3Ch))/nPartIm2;
GraphData.BarDataIm2(1,3) = sum((OutputData.distTableIm2(1,:) > OutputData.SAAcutoff3Ch) & (OutputData.distTableIm2(2,:) <= OutputData.SAAcutoff3Ch))/nPartIm2;
GraphData.BarDataIm2(1,4) = 1-sum(GraphData.BarDataIm2(1,1:3));

GraphData.BarDataIm2(2,1) = sum((rT2(1,:) <= OutputData.SAAcutoff3Ch) & (rT2(2,:) <= OutputData.SAAcutoff3Ch))/(Conditions.Iterations*nPartIm2);
GraphData.BarDataIm2(2,2) = sum((rT2(1,:) <= OutputData.SAAcutoff3Ch) & (rT2(2,:) > OutputData.SAAcutoff3Ch))/(Conditions.Iterations*nPartIm2);
GraphData.BarDataIm2(2,3) = sum((rT2(1,:) > OutputData.SAAcutoff3Ch) & (rT2(2,:) <= OutputData.SAAcutoff3Ch))/(Conditions.Iterations*nPartIm2);
GraphData.BarDataIm2(2,4) = 1-sum(GraphData.BarDataIm2(2,1:3));
clear rT2

rT3(1,:) = Randomizations.randTableIm1_2(:);
rT3(2,:) = Randomizations.randTableIm1_3(:);

GraphData.BarDataIm3(1,1) = sum((OutputData.distTableIm3(1,:) <= OutputData.SAAcutoff3Ch) & (OutputData.distTableIm3(2,:) <= OutputData.SAAcutoff3Ch))/nPartIm3;
GraphData.BarDataIm3(1,2) = sum((OutputData.distTableIm3(1,:) <= OutputData.SAAcutoff3Ch) & (OutputData.distTableIm3(2,:) > OutputData.SAAcutoff3Ch))/nPartIm3;
GraphData.BarDataIm3(1,3) = sum((OutputData.distTableIm3(1,:) > OutputData.SAAcutoff3Ch) & (OutputData.distTableIm3(2,:) <= OutputData.SAAcutoff3Ch))/nPartIm3;
GraphData.BarDataIm3(1,4) = 1-sum(GraphData.BarDataIm3(1,1:3));

GraphData.BarDataIm3(2,1) = sum((rT3(1,:) <= OutputData.SAAcutoff3Ch) & (rT3(2,:) <= OutputData.SAAcutoff3Ch))/(Conditions.Iterations*nPartIm3);
GraphData.BarDataIm3(2,2) = sum((rT3(1,:) <= OutputData.SAAcutoff3Ch) & (rT3(2,:) > OutputData.SAAcutoff3Ch))/(Conditions.Iterations*nPartIm3);
GraphData.BarDataIm3(2,3) = sum((rT3(1,:) > OutputData.SAAcutoff3Ch) & (rT3(2,:) <= OutputData.SAAcutoff3Ch))/(Conditions.Iterations*nPartIm3);
GraphData.BarDataIm3(2,4) = 1-sum(GraphData.BarDataIm3(2,1:3));
clear rT3


% Bar Graphs
h = figure(2);
tPos = get(figure(2), 'Position');
set(figure(2), 'Position', [75, tPos(2), 1120, 1450]);
set (figure(2), 'Color', 'w');
subplot (3,3,1);
bar(GraphData.BarDataIm1, 1, 'grouped');
colormap ('Gray');
legend ([Conditions.Ch1Label ':' Conditions.Ch2Label ':' Conditions.Ch3Label], [Conditions.Ch1Label ':' Conditions.Ch2Label], [Conditions.Ch1Label ':' Conditions.Ch3Label], 'Not Associated', 'Location', 'Best');
ylabel ('Fraction of Colocalized Particles');
title (['Colocalization of ' Conditions.Ch1Label]);
set(gca,'XTickLabel',{'Data', 'Siml. Random Positions'});

subplot (3,3,2);
bar(GraphData.BarDataIm2, 1, 'grouped');
colormap ('Gray');
legend ([Conditions.Ch2Label ':' Conditions.Ch1Label ':' Conditions.Ch3Label], [Conditions.Ch2Label ':' Conditions.Ch1Label], [Conditions.Ch2Label ':' Conditions.Ch3Label], 'Not Associated', 'Location', 'Best');
ylabel ('Fraction of Colocalized Particles');
title (['Colocalization of ' Conditions.Ch2Label]);
set(gca,'XTickLabel',{'Data', 'Siml. Random Positions'});

subplot (3,3,3);
bar(GraphData.BarDataIm3, 1, 'grouped');
colormap ('Gray');
legend ([Conditions.Ch3Label ':' Conditions.Ch1Label ':' Conditions.Ch2Label], [Conditions.Ch3Label ':' Conditions.Ch1Label], [Conditions.Ch3Label ':' Conditions.Ch2Label], 'Not Associated', 'Location', 'Best');
ylabel ('Fraction of Colocalized Particles');
title (['Colocalization of ' Conditions.Ch3Label]);
set(gca,'XTickLabel',{'Data', 'Siml. Random Positions'});
text(-2,1.2,['SAA Quantification: ' Conditions.Ch1Label ', ' Conditions.Ch2Label ' & ' Conditions.Ch3Label]);
saveas (figure(2), 'SAA - bargraph');
save2pdf ('SAA - bargraph.pdf', figure(2), 600);
close (h);

%Density Graphs
h = figure(3); 
set (figure(3), 'Color', 'w');
tPos = get(figure(3), 'Position');
set(figure(3), 'Position', [75, tPos(2), 1120, 1450]);
subplot (3,2,1);
contour(GraphData.Im1ContourFiltered, 100);
colormap ('Jet');
axis ('square');
xlabel (['Separation - ' Conditions.Ch1Label ':' Conditions.Ch2Label '(nm)']);
ylabel (['Separation - ' Conditions.Ch1Label ':' Conditions.Ch3Label '(nm)']);
title (Conditions.Ch1Label);

subplot (3,2,2);
contour(Randomizations.RandGraph1ContouredFiltered, 100);
colorbar('location','eastoutside');
axis ('square');
xlabel (['Separation - ' Conditions.Ch1Label ':' Conditions.Ch2Label '(nm)']);
ylabel (['Separation - ' Conditions.Ch1Label ':' Conditions.Ch3Label '(nm)']);
title ([Conditions.Ch1Label ' - Randomized']);

subplot (3,2,3);
contour(GraphData.Im2ContourFiltered, 100);
axis ('square');
xlabel (['Separation - ' Conditions.Ch2Label ':' Conditions.Ch1Label '(nm)']);
ylabel (['Separation - ' Conditions.Ch2Label ':' Conditions.Ch3Label '(nm)']);
title (Conditions.Ch2Label);

subplot (3,2,4);
contour(Randomizations.RandGraph2ContouredFiltered, 100);
colorbar('location','eastoutside')
axis ('square');
xlabel (['Separation - ' Conditions.Ch2Label ':' Conditions.Ch1Label '(nm)']);
ylabel (['Separation - ' Conditions.Ch2Label ':' Conditions.Ch3Label '(nm)']);
title ([Conditions.Ch2Label ' - Randomized']);

subplot (3,2,5);
contour(GraphData.Im3ContourFiltered, 100);
axis ('square');
xlabel (['Separation - ' Conditions.Ch3Label ':' Conditions.Ch1Label '(nm)']);
ylabel (['Separation - ' Conditions.Ch3Label ':' Conditions.Ch2Label '(nm)']);
title (Conditions.Ch3Label);

subplot (3,2,6);
contour(Randomizations.RandGraph3ContouredFiltered, 100);
colorbar('location','eastoutside')
axis ('square');
xlabel (['Separation - ' Conditions.Ch3Label ':' Conditions.Ch1Label '(nm)']);
ylabel (['Separation - ' Conditions.Ch3Label ':' Conditions.Ch2Label '(nm)']);
title ([Conditions.Ch3Label ' - Randomized']);

%add sub-plots
new_axis = axes('position',[ 0.31 0.825 .1 .1]); axis ('square'); axis([0,100,0,100]);
contour(GraphData.Im1ContourFiltered, 25);
axis ([0 100 0 100]);
axis ('square')
line ([OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], [0, Conditions.BinMax], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
line ([0, Conditions.BinMax], [OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');

new_axis = axes('position',[ 0.725 0.825 .1 .1]); axis ('square'); axis([0,100,0,100]);
contour(Randomizations.RandGraph1ContouredFiltered, 25);
axis ([0 100 0 100]);
axis ('square')
line ([OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], [0, Conditions.BinMax], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
line ([0, Conditions.BinMax], [OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');

new_axis = axes('position',[ 0.31 0.525 .1 .1]); axis ('square'); axis([0,100,0,100]);
contour(GraphData.Im2ContourFiltered, 25);
axis ([0 100 0 100]);
axis ('square')
line ([OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], [0, Conditions.BinMax], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
line ([0, Conditions.BinMax], [OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');

new_axis = axes('position',[ 0.725 0.525 .1 .1]); axis ('square'); axis([0,100,0,100]);
contour(Randomizations.RandGraph2ContouredFiltered, 25);
axis ([0 100 0 100]);
axis ('square')
line ([OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], [0, Conditions.BinMax], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
line ([0, Conditions.BinMax], [OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');

new_axis = axes('position',[ 0.31 0.225 .1 .1]); axis ('square'); axis([0,100,0,100]);
contour(GraphData.Im3ContourFiltered, 25);
axis ([0 100 0 100]);
axis ('square')
line ([OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], [0, Conditions.BinMax], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
line ([0, Conditions.BinMax], [OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');

new_axis = axes('position',[ 0.725 0.225 .1 .1]); axis ('square'); axis([0,100,0,100]);
contour(Randomizations.RandGraph3ContouredFiltered, 25);
axis ([0 100 0 100]);
axis ('square')
line ([OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], [0, Conditions.BinMax], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
line ([0, Conditions.BinMax], [OutputData.SAAcutoff3Ch, OutputData.SAAcutoff3Ch], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
hold off;
text(-350,880,['Three Channel SAA Density Maps: ' Conditions.Ch1Label ', ' Conditions.Ch2Label ' & ' Conditions.Ch3Label]);
saveas (figure(3), 'SAA - Density Plots');
save2pdf ('SAA - Density Plots.pdf', figure(3), 600);
close (h);

%Export data to gui
SAAstruct  = OutputData;

%% Cleanup
if openPool == 1 %close pool of opened by function
    delete (gcp);
end

display (' ');
display (['   ...SAA Processing complete in ' num2str(toc(ttime)) ' seconds.']);

end %end function

%% Sub-functions

% function to find nearest neighbour to pXY in x/y/z coordinate list PosLlist
function [n1, n2] = nearestNeighbour(pXY, PosList1, PosList2, Range)
[~, tCol1] = find(PosList1(1,:)>=(pXY(1)-Range) & PosList1(1,:) <= (pXY(1)+Range) & PosList1(2,:)>=(pXY(2)-Range) & PosList1(2,:)<=(pXY(2)+Range));
[~, tCol2] = find(PosList2(1,:)>=(pXY(1)-Range) & PosList2(1,:) <= (pXY(1)+Range) & PosList2(2,:)>=(pXY(2)-Range) & PosList2(2,:)<=(pXY(2)+Range));
if tCol1
    n1 = min(sqrt(bsxfun(@plus,dot(pXY,pXY,1)',dot(PosList1,PosList1,1))-2*pXY'*PosList1));
else
    n1 = NaN;
end

if tCol2
    n2 = min(sqrt(bsxfun(@plus,dot(pXY,pXY,1)',dot(PosList2,PosList2,1))-2*pXY'*PosList2));
else
    n2 = NaN;
end
end

function save2pdf(pdfFileName,handle,dpi)

% Backup previous settings
prePaperType = get(handle,'PaperType');
prePaperUnits = get(handle,'PaperUnits');
preUnits = get(handle,'Units');
prePaperPosition = get(handle,'PaperPosition');
prePaperSize = get(handle,'PaperSize');

% Make changing paper type possible
set(handle,'PaperType','<custom>');

% Set units to all be the same
set(handle,'PaperUnits','inches');
set(handle,'Units','inches');

% Set the page size and position to match the figure's dimensions
paperPosition = get(handle,'PaperPosition');
position = get(handle,'Position');
set(handle,'PaperPosition',[0,0,position(3:4)]);
set(handle,'PaperSize',position(3:4));

% Save the pdf (this is the same method used by "saveas")
print(handle,'-dpdf',pdfFileName,sprintf('-r%d',dpi))

% Restore the previous settings
set(handle,'PaperType',prePaperType);
set(handle,'PaperUnits',prePaperUnits);
set(handle,'Units',preUnits);
set(handle,'PaperPosition',prePaperPosition);
set(handle,'PaperSize',prePaperSize);
end