function [rvalue, g] = grFun( x,A,Start,End,Step,size_ROI)
    %Ripley Summary of this function goes here
    %   This Function calculates the ripley value for the data set x.
    %   x : nx2 matrix, n points, column 1 : x, column 2 : y
    %   A area of the ROI of interest
    %   Start : first search raduis
    %   End : last search raduis
    %   Step : step of search raduis
    %     x=1024*rand(1,10000);
    %     y=1024*rand(1,10000);

    N = size(x, 1);
    xmin = min(x, [], 1);
    xcenter = x(:,1)-xmin(1);
    ycenter = x(:,2)-xmin(2);
    i_range = (End - Start) / Step;
    i_range = i_range + 1;
    
    % If there's a parallel pool open, use that
    % Otherwise start one
    if isunix
            c = gcp('nocreate');
            if isempty(c)
                tempdir = getenv('TMPDIR');
                c = parcluster;
                c.JobStorageLocation = tempdir;
                c.NumWorkers=20;
                %parpool(c,4);
                parpool(c);
            end
    end
    
    poolobj = gcp;
    if ~any(~cellfun(@isempty, (strfind(poolobj.AttachedFiles(), 'compute_circle_area.m'))))
        addAttachedFiles(poolobj,{['ClusDoCUtils' filesep 'compute_circle_area.m']})
    end
    
    Kr = zeros(i_range,1);
    
    % for i=1:i_range 
    parfor i=1:i_range

        %i=i+1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute next step in radius
        r = (i-1)*Step;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        relative_area = compute_circle_area(xcenter',ycenter',r,size_ROI);
        %indexborder=find(relative_area<1);

        [idx, ~] = rangesearch(x,x,r); % can this be sped up?
        % Pontentially with pdist + binary
        % operations?

        Kfuncans1 = (cellfun('length',idx)-1)./relative_area;     % remove the identity
        Kr(i) = mean((Kfuncans1*A/N));    % ? 
        rvalue(i,1) = r;
    end
    %
    g = diff(Kr);
    g = [g; g(length(g))];
    g = g/Step./(2*pi*rvalue);
    %
end
