%==========================================================================
%                              FUNCTION
% Performs a hierarchal segmentation on the output of OPTICS.m in order to
% identify clusters, based on the RD plot.
%
%==========================================================================
%                             INPUT/OUTPUT
% Inputs:
%   -InputOPTICS: OPTICStable matrix produced by OPTICS.m
%   -minClust: minimum number of molecules to be considered a cluster
%   -minSplit: minimum difference between a peak used for a split and the
%              mean RD value of the resulting segments for those segments
%              to be considered a cluster. A value of 0.75 is default.
%
% Output:OPTICS ANALYSIS
%   -Cluster(ii): a structure containing ii clusters, where:
%       -.parentMat: RD table for parental cluster
%       -.l/rCluster: RD table for the left/right cluster of the current 
%                     split
%       -.l/rRatio: peak-to-mean RD ratio of the left/right segment in the 
%                 current split
%       -.Parent: number of the cluster which is the parent to the current 
%                 split
%       -.r/lProcessed: Whether the left/right segment is further 
%                       processed: 1 = split further, 3 = leaf
%       -.l/rChild: number of the cluster for the left/right child of the
%                   current split
%       -.l/rEdgePts: for 2D images, provides x/y coordinates of the
%                     bounding edge of the cluster; for 3D images provides
%                     vertices information for use with trisurf
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

function Cluster = hierOPTICS(InputOPTICS,minClust,minSplit)

%% Check inputs
if (nargin ~= 2) && (nargin ~= 3)
    edit ('hierOPTICS.m');
    error ('hierOPTICS must be called with 2 or 3 input variables');
end

if size (InputOPTICS,2) ~= 6
    error ('Input must be the 6-column ''OPTICStable'' matrix produced by OPTICS.m');
end

if minClust < 2
    display ('Minimum cluster size is too small, setting to 2.');
    minClust = 2;
end

if (minSplit <= 0) || (minSplit >= 1)
    error ('Minsplit needs to be between 0 and 1');
end

if nargin == 2
    minSplit = 0.75;
end

isZ = max(InputOPTICS(:,6)); %check to see if there is z-dimention

%% Segment RD plot
% Set first cluster as the root
[lCluster, lRatio, rCluster, rRatio] = splitOPTICS (InputOPTICS);

Cluster(1).parentMat = InputOPTICS;
Cluster(1).lCluster = lCluster;
Cluster(1).lRatio = lRatio;
Cluster(1).rCluster = rCluster;
Cluster(1).rRatio = rRatio;
Cluster(1).Parent = 0; %is root
Cluster(1).rProcessed = 0; %flag - right child processed
Cluster(1).lProcessed = 0; %flag - left child processed
Cluster(1).lChild = 1;
Cluster(1).rChild = 2;

Done = 0; % flag
nClusters = 1; %number of split Clusters

while ~Done
   %find next Cluster to process
   pClust = 1;  %marker
   
   %process all existing Clusters 
   while (pClust <= nClusters)
        if ~Cluster(pClust).lProcessed %if left Cluster is unprocessed
            [lCluster, lRatio, rCluster, rRatio] = splitOPTICS (Cluster(pClust).lCluster);
            if (size(lCluster,1) >= minClust) && (size(rCluster,1) >= minClust) && (lRatio <= minSplit) && (rRatio <= minSplit) %is real split
                Cluster(pClust).lChild = nClusters + 1;
                Cluster(pClust).lProcessed = 1;
                Cluster(nClusters+1).parentMat = Cluster(pClust).lCluster;
                Cluster(nClusters+1).lCluster = lCluster;
                Cluster(nClusters+1).lRatio = lRatio;
                Cluster(nClusters+1).rCluster = rCluster;
                Cluster(nClusters+1).rRatio = rRatio;
                Cluster(nClusters+1).Parent = pClust; %parent
                Cluster(nClusters+1).rProcessed = 0; %flag - right child processed
                Cluster(nClusters+1).lProcessed = 0; %flag - left child processed
                nClusters = nClusters + 1;

                %find boundary points
                try
                    if isZ %if 3D image
                        Cluster(pClust).lEdgePts = convhull (Cluster(pClust).lCluster(:,4:6));
                    else %is 2D image
                        edgePts = convhull (Cluster(pClust).lCluster(:,4:5));
                        Cluster(pClust).lEdgePts = Cluster(pClust).lCluster(edgePts,4:5);
                    end
                catch
                    Cluster(pClust).lEdgePts = NaN;
                end
                
            elseif size(lCluster,1) >= minClust && (lRatio <= minSplit) %left Cluster is sufficiently big
                Cluster(pClust).lCluster = lCluster; %re-write current Cluster back to parent
                Cluster(pClust).lProcessed = 0;
                
            elseif size(rCluster,1) >= minClust && (rRatio <= minSplit) %right Cluster is sufficiently big
                Cluster(pClust).lCluster = rCluster; %re-write current Cluster back to parent
                Cluster(pClust).lProcessed = 0;

            else %parent is leaf
                Cluster(pClust).lProcessed = 3;
            end
        end
            
            %Process right Cluster
        if ~Cluster(pClust).rProcessed %if right Cluster is unprocessed
            [lCluster, lRatio, rCluster, rRatio] = splitOPTICS (Cluster(pClust).rCluster);
            if (size(lCluster,1) >= minClust) && (size(rCluster,1) >= minClust) && (lRatio <= minSplit) && (rRatio <= minSplit) %is real split
                Cluster(pClust).rChild = nClusters + 1;
                Cluster(pClust).rProcessed = 1;
                Cluster(nClusters+1).parentMat = Cluster(pClust).rCluster;
                Cluster(nClusters+1).lCluster = lCluster;
                Cluster(nClusters+1).lRatio = lRatio;
                Cluster(nClusters+1).rCluster = rCluster;
                Cluster(nClusters+1).rRatio = rRatio;
                Cluster(nClusters+1).Parent = pClust; %parent
                Cluster(nClusters+1).rProcessed = 0; %flag - right child processed
                Cluster(nClusters+1).lProcessed = 0; %flag - left child processed
                nClusters = nClusters + 1;

                %find boundary points
                try
                    if isZ %if 3D image
                        Cluster(pClust).rEdgePts = convhull (Cluster(pClust).rCluster(:,4:6));
                    else %is 2D image
                        edgePts = convhull (Cluster(pClust).rCluster(:,4:5));
                        Cluster(pClust).rEdgePts = Cluster(pClust).rCluster(edgePts,4:5);
                    end
                catch
                    Cluster(pClust).lEdgePts = NaN;
                end
                
            elseif size(lCluster,1) >= minClust && (lRatio <= minSplit) %left Cluster is sufficiently big
                Cluster(pClust).rCluster = lCluster; %re-write current Cluster back to parent
                Cluster(pClust).rProcessed = 0;
                
            elseif size(rCluster,1) >= minClust && (rRatio <= minSplit) %right Cluster is sufficiently big
                Cluster(pClust).rCluster = rCluster; %re-write current Cluster back to parent
                Cluster(pClust).rProcessed = 0;

            else %parent is leaf
                Cluster(pClust).rProcessed = 3;
            end
        end
        pClust = pClust + 1;
    end %end current Cluster processig
   
    Done = 1;
    for ii=1:length(Cluster)
        if (Cluster(ii).lProcessed == 0) || (Cluster(ii).rProcessed == 0)
            Done = 0;
        end
    end

end

end



function [lCluster, lRatio, rCluster, rRatio] = splitOPTICS (InputOPTICS)
% function to split an OPTICS table (from OPTICS.m) into two Clusters;
% centred on the highest RD peak in the dataset.

%Inputs: 
%   -InputOPTICS = OPTICS table from OPTICS.m

%Outputs:
%   -l/rCluster: OPTICS table formatted left/right Clusters
%   l/rRatio: Ratio of mean RD value for splits vs. RD value of split
%   point

Split = find(InputOPTICS(:,2) == max(InputOPTICS(:,2)),1);
SplitV = InputOPTICS(Split,2);
lCluster = InputOPTICS(1:(Split-1),:);
lRatio = mean(lCluster(:,2))/SplitV;
rCluster = InputOPTICS((Split+1):end,:);
rRatio = mean(rCluster(:,2))/SplitV;

end