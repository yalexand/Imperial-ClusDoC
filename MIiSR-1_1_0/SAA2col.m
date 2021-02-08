function SAAstruct = SAA2col(CroppedPos, nCores, Conditions)

%==========================================================================
%                              FUNCTION
% Function to perform spatial association analysis on GSD position files
% using a position/intensity file produced from GSDLoadCrop.
%
% Function which performs a spatial association analysis of
% single-molecules images.  The Function measures the distances between the 
% particles in Im1 and the closest particle in Im2, and then compares 
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
%             cutoff, p = 2.00 for 95% cut-off
%       -Cutoff = value to add to the cut-off determined by function.  Used
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

    %SAA Calculations
    Conditions.AnalyzedFraction = 1; %fraction of dataset to analyze
    Conditions.SD = 1.65; %Standard deviation of error.  1.65 = p of 0.90; 2.00 = p of 0.95.
    Conditions.Cutoff = 0; %Value to add to SAA, for known separations or poor image registration
    Conditions.BinMax = 1000; %maximum distance to assess
    Conditions.Iterations = 5; %randomization iterations for Simulated Random Positions calculations. Minimum of 1.
    Conditions.polyRand = 0; %Perform randomization over an area deterined using convhull; used to ensure that clear 
                             % areas created by invisible objects are not
                             % included in SAA calculations. 0 = off; 1 = on
end

% Outputs:
%
%   -OutputData
%       .Im1MeaPhotons/Im2: Average photons in Im1/2
%       .Im1MeanPrecission/Im2: Mean precision of Im1/2.
%           -Calculated as 200/sqrt(number photons)
%       .sigmaRMS: Root-mean squared of Im1/Im2 precision
%       .SAAcutoff: cutoff for SAA analysis, 
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
%       .Im1/Im2 - Image 1/2 file name
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
    open ('GSD_SAA2d.m');
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

if Conditions.Iterations < 1
    display ('Conditions.Iterations must be 1 or greater. Setting Conditions.Iterations to 1.');
    Conditions.Iterations = 1;
end

%% Load Data & Prepare for Processing
RawData.Conditions = Conditions;
xAxisLabels = [0:Conditions.BinMax/100:Conditions.BinMax];

load (CroppedPos);
if exist ('cropPos', 'var')
    RawData.Im1Pos = cropPos.Ch(1).Pos;
    RawData.Im2Pos = cropPos.Ch(2).Pos;
    clear cropPos
else
    RawData.Im1Pos = allPos.Ch(1).Pos;
    RawData.Im2Pos = allPos.Ch(2).Pos;
    clear allPos
end

%reduce data to subset, if required
if Conditions.AnalyzedFraction < 1
    aT = tic;
    display (' ');
    display ('Reducing data size.');
    PosKeep = rand(length(RawData.Im1Pos),1);
    RawData.Im1Pos = RawData.Im1Pos(PosKeep <= Conditions.AnalyzedFraction,:);
    PosKeep = rand(length(RawData.Im2Pos),1);
    RawData.Im2Pos = RawData.Im2Pos(PosKeep <= Conditions.AnalyzedFraction,:);
    display (['. . .Data reduction completed in ' num2str(toc(aT)) ' seconds.']);
    clear PosKeep aT;
end %end subset generation

%determine if data is 2D or 3D
nDims = max([max(RawData.Im1Pos(:,3)), max(RawData.Im2Pos(:,3))]);
if nDims ~= 0
    nDims = 3;
else
    nDims = 2;
end

%Calculate precision & sigma-RMS
OutputData.Im1MeanPhotons = mean(RawData.Im1Pos(:,4));
OutputData.Im1MeanPrecission = mean(RawData.Im1Pos(:,5));
OutputData.Im2MeanPhotons = mean(RawData.Im2Pos(:,4));
OutputData.Im2MeanPrecission = mean(RawData.Im2Pos(:,5));
OutputData.sigmaRMS = sqrt((OutputData.Im1MeanPrecission^2) + (OutputData.Im2MeanPrecission^2));
OutputData.SAAcutoff = (OutputData.sigmaRMS*Conditions.SD) + Conditions.Cutoff;

%% Generate Distance Tables
% Using fastest eucld. distance equation - http://www.mathworks.com/matlabcentral/fileexchange/71-distance-m
display ('...Performing 2-Channel SAA Analysis ');
display ('   ...Generating Distance Tables.');
nPartIm1 = size(RawData.Im1Pos,1);
nPartIm2 = size(RawData.Im2Pos,1);
Range = Conditions.BinMax;

%pre-allocate for speed
distTableIm1(nPartIm1,1) = double(0);
distTableIm2(nPartIm2,1) = double(0);
Pos1 = RawData.Im1Pos(:,1:3)';
Pos2 = RawData.Im2Pos(:,1:3)';

aT = tic;
parfor f=1:nPartIm1
    distTableIm1(f) =  nearestNeighbour(Pos1(:,f), Pos2, Range);
end
%Data for Graphs
OutputData.distTableIm1 = distTableIm1';
OutputData.Im1Graph = hist(distTableIm1', xAxisLabels);
OutputData.Im1Graph = OutputData.Im1Graph./nPartIm1;
clear distTableIm1;
display ('      -Distance table calculated for dataset 1.');
display (['          -completed in ' num2str(toc(aT)) ' seconds.']);
clear aT

aT = tic;
parfor f=1:nPartIm2
    distTableIm2(f) =  nearestNeighbour(Pos2(:,f), Pos1, Range);
end
%Data for Graphs
OutputData.distTableIm2 = distTableIm2';
OutputData.Im2Graph = hist(distTableIm2', xAxisLabels);
OutputData.Im2Graph = OutputData.Im2Graph./nPartIm2;
clear distTableIm2;
display ('      -Distance table calculated for dataset 2.');
display (['          -completed in ' num2str(toc(aT)) ' seconds.']);
clear aT

%% Random Populations
if Conditions.Iterations
    display (' ');
    display ('   ...Generating Randomized Distance Tables.');
    
    aT = tic;

    % determine image area to use for randomization
    if Conditions.polyRand % if ROI is refined using convenx hull
        if nDims == 2
            tPos = vertcat(RawData.Im1Pos(:,1:2), RawData.Im2Pos(:,1:2));
            [~, tPos] = convhull(tPos); %area of convex hull
            dX = sqrt(tPos);
            dY = dX;
            dZ = 0;
            clear tPos
        else % if using user-set rectangular area
            tPos = vertcat(RawData.Im1Pos(:,1:3), RawData.Im2Pos(:,1:3));
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
    end
    clear dX dY dZ RandSize*;
    
    
    for i=1:Conditions.Iterations
        tTable2 = randTable2(i).a;
        parfor f=1:nPartIm1
            randTableIm1(i,f) =  nearestNeighbour(Pos1(:,f), tTable2, Range);
        end
    end
    %data for graphs
    OutputData.randTableIm1 = randTableIm1';
    OutputData.Im1RandGraph = hist(OutputData.randTableIm1, xAxisLabels);
    for i=1:Conditions.Iterations
        OutputData.Im1RandGraph(:,i) = OutputData.Im1RandGraph(:,i)./nPartIm1;
    end
    OutputData.Im1RandGraph(:,Conditions.Iterations+1) = mean(OutputData.Im1RandGraph, 2);
    clear randTableIm1;
    display (['      -' num2str(Conditions.Iterations) ' randomized distance tables calculated for dataset 1.']);
    display (['          -completed in ' num2str(toc(aT)) ' seconds.']);
    clear aT 

    aT = tic;   
    for i=1:Conditions.Iterations
        tTable1 = randTable1(i).a;
        parfor f=1:nPartIm2
            randTableIm2(i,f) =  nearestNeighbour(Pos2(:,f), tTable1, Range);
        end
    end
    %data for graphs
    OutputData.randTableIm2 = randTableIm2';
    OutputData.Im2RandGraph = hist(OutputData.randTableIm2, xAxisLabels);
    for i=1:Conditions.Iterations
        OutputData.Im2RandGraph(:,i) = OutputData.Im2RandGraph(:,i)./nPartIm2;
    end
    OutputData.Im2RandGraph(:,Conditions.Iterations+1) = mean(OutputData.Im2RandGraph, 2);
    display (['      -' num2str(Conditions.Iterations) ' randomized distance tables calculated for dataset 2.']);
    display (['          -completed in ' num2str(toc(aT)) ' seconds.']);
    clear aT tTable1 tTable2 randTableIm2 randTable1 randTable2;

end

clear Pos1 Pos2; 

%% Graphs
OutputData.xAxisLabels = xAxisLabels;

%Plot Im1:Im2
Im1Max = max(OutputData.Im1Graph);
Im1Max(2) = max(OutputData.Im1RandGraph(:,end));
Im1Max = max(ceil(Im1Max.*20)./20); %round up to nearest 0.05

h = figure(1);
plot (OutputData.xAxisLabels, OutputData.Im1Graph,'k');
hold on;
plot (OutputData.xAxisLabels, OutputData.Im1RandGraph(:,end), '--k');
line ([OutputData.SAAcutoff, OutputData.SAAcutoff], [0, Im1Max], 'LineStyle', '-.', 'Color', 'k');
xlabel ([Conditions.Ch1Label ':' Conditions.Ch2Label ' Separation Distance']);
ylabel ('Fraction of Particles');
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
axis ([0, Conditions.BinMax, 0, Im1Max]);
set (figure(1), 'Color', 'w');
hold off;
saveas (figure(1), ['SAA - ' Conditions.fName ' - ' Conditions.Ch1Label ' vs ' Conditions.Ch2Label '.fig']);
pause(0.5);
close (h);

%Plot Im2:Im1
Im2Max = max(OutputData.Im2Graph);
Im2Max(2) = max(OutputData.Im2RandGraph(:,end));
Im2Max = max(ceil(Im2Max.*20)./20); %round up to nearest 0.05

h = figure(2);
plot (OutputData.xAxisLabels, OutputData.Im2Graph,'k');
hold on;
plot (OutputData.xAxisLabels, OutputData.Im2RandGraph(:,end), '--k');
line ([OutputData.SAAcutoff, OutputData.SAAcutoff], [0, Im2Max], 'LineStyle', '-.', 'Color', 'k');
xlabel ([Conditions.Ch2Label ':' Conditions.Ch1Label ' Separation Distance']);
ylabel ('Fraction of Particles');
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
axis ([0, Conditions.BinMax, 0, Im1Max]);
set (figure(2), 'Color', 'w');
hold off;
saveas (figure(2), ['SAA - ' Conditions.fName ' - '  Conditions.Ch2Label ' vs ' Conditions.Ch1Label '.fig']);
pause(0.5);
close (h);

%Plot Bar Data
OutputData.BarData(1,1) = sum(OutputData.distTableIm1 <= OutputData.SAAcutoff);
OutputData.BarData(1,2) = sum(OutputData.randTableIm1(:) <= OutputData.SAAcutoff)/Conditions.Iterations;
OutputData.BarData(2,1) = sum(OutputData.distTableIm2 <= OutputData.SAAcutoff);
OutputData.BarData(2,2) = sum(OutputData.randTableIm2(:) <= OutputData.SAAcutoff)/Conditions.Iterations;

OutputData.eBarData(1) = OutputData.BarData(1,1)./OutputData.BarData(1,2);
OutputData.eBarData(2) = OutputData.BarData(2,1)./OutputData.BarData(2,2);
OutputData.fBarData(1,1) = (sum(OutputData.distTableIm1 <= OutputData.SAAcutoff))/nPartIm1;
OutputData.fBarData(1,2) = (sum(OutputData.randTableIm1(:) <= OutputData.SAAcutoff)/Conditions.Iterations)/nPartIm1;
OutputData.fBarData(2,1) = (sum(OutputData.distTableIm2 <= OutputData.SAAcutoff))/nPartIm2;
OutputData.fBarData(2,2) = (sum(OutputData.randTableIm2(:) <= OutputData.SAAcutoff)/Conditions.Iterations)/nPartIm2;

h = figure(3);
set (figure(3), 'Color', 'w');
subplot (1,2,1);
bar(OutputData.BarData, 1, 'grouped');
colormap ('Gray');
legend ('Data', 'Sim. Positions');
ylabel ('Number of Colocalized Particles');
set(gca,'XTickLabel',{[Conditions.Ch1Label ':' Conditions.Ch2Label];[Conditions.Ch2Label ':' Conditions.Ch1Label]});
subplot (1,2,2);
bar(OutputData.fBarData, 1, 'grouped');
colormap ('Gray');
ylabel ('Fraction of Colocalized Particles');
set(gca,'XTickLabel',{[Conditions.Ch1Label ':' Conditions.Ch2Label];[Conditions.Ch2Label ':' Conditions.Ch1Label]});
saveas (figure(3), ['SAA - ' Conditions.fName ' - Fraction Interacting of ' Conditions.Ch1Label ' and ' Conditions.Ch2Label '.fig']);
pause(0.5);
close (h);

h = figure(4);
bar([0,2], OutputData.eBarData, 'k');
ylabel ('Fold-Colocalization (Rel. Sim. Positions)');
set(gca,'XTickLabel',{[Conditions.Ch1Label ':' Conditions.Ch2Label];[Conditions.Ch2Label ':' Conditions.Ch1Label]});
set (figure(4), 'Color', 'w');
saveas (figure(4), ['SAA - ' Conditions.fName ' - Fold Enrichment of ' Conditions.Ch1Label ' and ' Conditions.Ch2Label '.fig']);
pause(0.5);
close (h);

%lab book figure
h = figure(5);
set (figure(5), 'Color', 'w');
subplot(2,2,1);
plot (OutputData.xAxisLabels, OutputData.Im1Graph,'k');
hold on;
plot (OutputData.xAxisLabels, OutputData.Im1RandGraph(:,end), '--k');
line ([OutputData.SAAcutoff, OutputData.SAAcutoff], [0, Im1Max], 'LineStyle', '-.', 'Color', 'k');
axis ([0, Conditions.BinMax, 0, Im1Max]);
xlabel ([Conditions.Ch1Label ':' Conditions.Ch2Label ' Separation Distance']);
ylabel ('Fraction of Particles');
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
hold off;

subplot(2,2,2);
plot (OutputData.xAxisLabels, OutputData.Im2Graph,'k');
hold on;
plot (OutputData.xAxisLabels, OutputData.Im2RandGraph(:,end), '--k');
line ([OutputData.SAAcutoff, OutputData.SAAcutoff], [0, Im2Max], 'LineStyle', '-.', 'Color', 'k');
axis ([0, Conditions.BinMax, 0, Im1Max]);
xlabel ([Conditions.Ch2Label ':' Conditions.Ch1Label ' Separation Distance']);
ylabel ('Fraction of Particles');
legend ('Data', 'Sim. Positions', 'SAA Cutoff');
axis ([0, Conditions.BinMax, 0, Im1Max]);
hold off;

subplot(2,2,3);
bar(OutputData.fBarData, 1, 'grouped');
colormap ('Gray');
legend ('Data', 'Sim. Positions');
ylabel ('Fraction of Colocalized Particles');
set(gca,'XTickLabel',{[Conditions.Ch1Label ':' Conditions.Ch2Label];[Conditions.Ch2Label ':' Conditions.Ch1Label]});

subplot(2,2,4);
bar([0,2], OutputData.eBarData, 'k');
ylabel ('Fold-Colocalization (Rel. Sim. Positions)');
set(gca,'XTickLabel',{[Conditions.Ch1Label ':' Conditions.Ch2Label];[Conditions.Ch2Label ':' Conditions.Ch1Label]});
tPos = get(figure(5), 'Position');
set(figure(5), 'Position', [(tPos(1)/2), (tPos(2)/2), 1120, 840]);
pause (0.5); %pause to allow for resize
saveas (figure(5), ['SAA Net Data - ' Conditions.fName ' - ' Conditions.Ch1Label ' vs ' Conditions.Ch2Label '.fig']);
saveas (figure(5), ['SAA Net Data - ' Conditions.fName ' - ' Conditions.Ch1Label ' vs ' Conditions.Ch2Label '.pdf']);
pause(0.5);
close (h);

%Send data to GUI
SAAstruct  = OutputData;

%% Clean-up
if openPool == 1 %close pool of opened by function
    delete (gcp);
end

display (' ');
display (['   ...SAA Processing complete in ' num2str(toc(ttime)) ' seconds.']);

end %end function

% function to find nearest neighbour to pXY in x/y/z coordinate list PosLlist
function nn = nearestNeighbour(pXY, PosList, Range)
[~, tCol] = find(PosList(1,:)>=(pXY(1)-Range) & PosList(1,:) <= (pXY(1)+Range) & PosList(2,:)>=(pXY(2)-Range) & PosList(2,:)<=(pXY(2)+Range));
if tCol
    nn = min(sqrt(bsxfun(@plus,dot(pXY,pXY,1)',dot(PosList,PosList,1))-2*pXY'*PosList));
else
    nn = NaN;
end
end