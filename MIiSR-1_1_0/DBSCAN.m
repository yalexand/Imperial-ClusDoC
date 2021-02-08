function [DBSCANmat, DBSCANtable, eps] = DBSCAN(CroppedPos, Ch1, k,Eps, plotClusters, fName)

%==========================================================================
%                              FUNCTION
% Performs DBSCAN analysis for cluster identification on molecule position
% file.
%
% Modified version of dbscan.m written by Michal Daszykowski,
% Department of Chemometrics, Institute of Chemistry, 
% The University of Silesia
% http://www.chemometria.us.edu.pl
%
%==========================================================================
%                             INPUT/OUTPUT
% Input: 
%   -CroppedPos: Position file produced by GSDLoadCrop
%   -Ch1: Which channel in CroppedPos to analyse
%   -k:minimum number of objects in a neighbourhood of an object (minimum
%      cluster size)
%   -Eps: neighbourhood radius, if not known avoid this parameter or put []
%   -plotClusters: generate an image with each cluster circled, 0 = no
%                  image, 1 = geneate image. Unclustered pixels in grey;
%                  clustered are colour-coded
%   -fName: name for saved images
%
% Output: 
%   -DBSCANmat: Structure containing:
%       -.k: k value entered by user
%       -.Eps: epsilon value used for analysis
%       -.Noise: x/y/z array of non-clustered (noise) molecules in the image
%       -.Cluster(ii): structure containing, where ii = cluster number:
%           -.Pos: x/y/z array of x/y/z positions of molecules in the
%                  cluster
%           -.edgePts: for 2D images, provides x/y coordinates of the
%                      bounding edge of the cluster; for 3D images provides
%                      vertice information for use with trisurf
%   -DBSCANtable - n x 5 table:
%       -Columns 1-3: x/y/z molecule positions
%       -Column 4: cluster ID number to which the molecule belongs. -1 =
%                  noise (unclustered) molecule
%       -Column 5: molecule type. core = 1, border = 0, noise/unclustered
%                  = -1)
%   -eps: Epsilon value used (either user-input 'Eps', or determined
%         automatically
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
%		Structures. 2015 PLoS Computational Biology
% 
% Please reference this paper in any publications which use this script for 
% analysis.
%==========================================================================
%% Check inputs & load file
if nargin < 3
    edit ('DBSCAN.m');
    error ('DBSCAN must be called with at least 3 input variables.');
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
    error ('Ch1 value is greater than the number of available channels');
end

%Extract x/y/z
x = allPos.Ch(Ch1).Pos(:,1:3);
clear allPos;

isZ = max(x(:,3));

%remove z coords for 2D data
if ~isZ
    x(:,3) = [];
end

x2 = x; %preserve original x-variable

[m,n] = size(x); %array size

%Define epsilon, if needed
if nargin == 3 || isempty(Eps)
    Eps=((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);
end

%% Perform DBSCAN analysis
eps = Eps;
x=[[1:m]' x];
type=zeros(1,m);
no=1;
touched=zeros(m,1);

for i=1:m
    if touched(i)==0;
       ob=x(i,:);
       D=dist2(ob(2:n),x(:,2:n));
       ind=find(D<=Eps);
    
       if length(ind)>1 && length(ind)<k+1       
          type(i)=0;
          class(i)=0;
       end
       if length(ind)==1
          type(i)=-1;
          class(i)=-1;  
          touched(i)=1;
       end

       if length(ind)>=k+1; 
          type(i)=1;
          class(ind)=ones(length(ind),1)*max(no);
          
          while ~isempty(ind)
                ob=x(ind(1),:);
                touched(ind(1))=1;
                ind(1)=[];
                D=dist2(ob(2:n),x(:,2:n));
                i1=find(D<=Eps);
     
                if length(i1)>1
                   class(i1)=no;
                   if length(i1)>=k+1;
                      type(ob(1))=1;
                   else
                      type(ob(1))=0;
                   end

                   for i=1:length(i1)
                       if touched(i1(i))==0
                          touched(i1(i))=1;
                          ind=[ind i1(i)];   
                          class(i1(i))=no;
                       end                    
                   end
                end
          end
          no=no+1; 
       end
   end
end

i1=find(class==0);
class(i1)=-1;
type(i1)=-1;

%% Generate Output
DBSCANmat.k = k;
DBSCANmat.Eps = eps;
DBSCANmat.Noise = x2(class==-1,:);
if max(class) == -1
    display ('No clusters identified.');
else
    for ii=1:max(class)
        DBSCANmat.Cluster(ii).Pos = x2(class==ii,:);
        try
            if n == 2 %2D
                edgePts = convhull (DBSCANmat.Cluster(ii).Pos);
                DBSCANmat.Cluster(ii).edgePts = DBSCANmat.Cluster(ii).Pos(edgePts,:);
            else %3D
                DBSCANmat.Cluster(ii).edgePts = convhull (DBSCANmat.Cluster(ii).Pos);
            end
        catch
            DBSCANmat.Cluster(ii).edgePts = NaN;
        end
    end
end

if n == 2
    x2(:,3) = 0;
end

DBSCANtable(:,1:3) = x2;
DBSCANtable(:,4) = class;
DBSCANtable(:,5) = type;

%% Plot image
if plotClusters && max(class) == -1
    h = figure(1);
    scatter (DBSCANmat.Noise(:,2), DBSCANmat.Noise(:,2), 1, [0.2 0.2 0.2]);
    axis ('equal'); 
    title (['DBSCAN, no clusters detected. k = ' num2str(k), ', ' char(949) ' = ' num2str(Eps)]);
    saveas (figure(1), ['DBSCAN - ' fName '.fig']);
    close (h);
elseif plotClusters
    clustercols = lines(max(class));
    h = figure(1);
    scatter (DBSCANmat.Noise(:,2), DBSCANmat.Noise(:,2), 1, [0.2 0.2 0.2]);
    hold on
    for ii=1:max(class)
        scatter(DBSCANmat.Cluster(ii).Pos(:,2), DBSCANmat.Cluster(ii).Pos(:,1), 1, clustercols(ii,:));
    end
    axis ('equal');
    title (['DBSCAN, ' num2str(max(class)) ' clusters detected. k = ' num2str(k), ', ' char(949) ' = ' num2str(Eps)]);
    hold off
    saveas (figure(1), ['DBSCAN - ' fName '.fig']);
    close (h);
end %end plotting

end


%% Sub-functions
function [D]=dist2(i,x)

% function: [D]=dist2(i,x)
%
% Calculates the Euclidean distances between the i-th object and all objects in x	 

[m,n]=size(x);
D=sqrt(sum((((ones(m,1)*i)-x).^2)'));

if n==1
   D=abs((ones(m,1)*i-x))';
end

end % end dist2