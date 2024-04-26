%-----------------------------------------------------------
function out = gr_assisted_quasi_cluster_analysis(obj, u, u0, s1, s2, A1, A2, save_dir_name)

params = obj.GR_ASSISTED_CLUSTERING;

% minimum #localizations in small and largr clusters
minPts1 = params.minPts1;
minPts2 = params.minPts2;

% if N3s/N2s (#locs within 3sigma divided by #locs within 2sigma) is way bigger than this number, it means that it is not standalone small cluster
t_ratio_3s2s = params.Small_Clusters_ratio_3s2s_threshold; 

    tic

    points = single(u);

    [sx,sy] = size(u);
    
    % BTW - coordinaets of all localizations
    % loc_coord = find(points == 1);
    
    %
    Ntot = sum(points(:));
    AREA = sum(u0,'All');
    rho = Ntot/AREA;
    % DEBUG
    
    % special proc shoul dbe written
    %
    %{
    When there are "no clusters" (only PSF), 
    or when large cluster size is unreasonably large, 
    so that it makes sense to disregard it, 
    or if one of the components comes out with too small amplitude, 
    or when 2 fitted sigmas are very close to each other - one needs to apply separate procedure parameterized by a single size.
    %}

    nloc_small = [];
    nloc_large = [];

    single_size_mode = (s2/s1) < params.min_Z2_Z1_ratio || s2 > params.max_Z2_allowed;

       if single_size_mode

           if (s2>s1) < 2
                s = (s1+s2)/2;
           elseif s2 > 500
                s = s1;
           end

           Np = gsderiv(points,s,0)*4*pi*s^2; % N multiplied by bell function of a height 1

           Np(Np<minPts1) = 0; % not interested in too sparse small clusters
           %  
           [maxima,maxima_v,maxima_ind,maxima_most_probable_v] = get_local_maxima(Np,sqrt(2)*s); 
           most_probable_clusters_N = maxima_most_probable_v;
           %            
           out = actual_point_clustering(points,maxima,[],s,s);
           %
           nloc_small = sample_nlocs(points, maxima, round(1.5*s));
 
       else

           norm1 = 1/(4*pi*s1^2);
           norm2 = 1/(4*pi*s2^2);
    
           Np1 = gsderiv(points,s1,0)/norm1; % N multiplied by bell function of a height 1
    
           Np1(Np1 < minPts1) = 0; % not interested in too sparse small clusters
           %  
           [maxima_small,maxima_small_v,maxima_small_ind,maxima_small_most_probable_v] = get_local_maxima(Np1,sqrt(2)*s1); % centres of small clusters
           most_probable_small_clusters_N = maxima_small_most_probable_v;
           %
           % peak vicinity analysis via 3-2 sigma collar counting
                   circle_2s1 = round_mask(2*s1);
                   circle_3s1 = round_mask(3*s1);
                   a3 = floor(length(circle_3s1)/2); % half-radius  
                   a2 = floor(length(circle_2s1)/2); % half-radius  
                   pad = a3 + 4; % specify the width of padding
                   padded_points = padarray(points, [pad, pad], 0, 'both');                    
                   x = maxima_small(:,1) + pad;
                   y = maxima_small(:,2) + pad;
                   N_2s = zeros(size(x)); 
                   N_3s = zeros(size(x));

                   for k=1:length(x)
                      xc=x(k);
                      yc=y(k);
                      rx3 = (xc-a3):(xc+a3);
                      ry3 = (yc-a3):(yc+a3); 
                      rx1 = (xc-a2):(xc+a2);
                      ry1 = (yc-a2):(yc+a2);  
                      %
                      u2 = padded_points(rx1,ry1).*circle_2s1;
                      u3 = padded_points(rx3,ry3).*circle_3s1;
                      N_2s(k) = sum(u2(:));
                      N_3s(k) = sum(u3(:));
                   end
                   % 
                   ratio_3s2s = N_3s./N_2s; % if there is no pedestal, this should be close to 1, otherwise bigger
                   %
                   exclusion_mask = ratio_3s2s > t_ratio_3s2s;
                   %
                   x_small_standalone = x(~exclusion_mask) - pad;
                   y_small_standalone = y(~exclusion_mask) - pad;

                   % corrected/redefined maxima_small
                   maxima_small = [x_small_standalone y_small_standalone];
                    
           % define large clusters' centres using density image from which contribution of small clusters was excluded
           where_small_clusters_are = imdilate(imdraw(sx,sy,maxima_small,1),strel('disk',round(2*s1)));
           Np2 = gsderiv(points & ~where_small_clusters_are,s2,0)/norm2;            
           Np2(Np2 < minPts2) = 0; % also don't take into account maxima that are too faint
           maxima_large = get_local_maxima(Np2,sqrt(2)*s2); 
           %
           % sampling
           Z1_sampling_fac = params.sampling_Z1_excess_factor;
           Z2_sampling_fac = params.sampling_Z2_excess_factor;
           %                  
           [out,perim] = actual_point_clustering(points,maxima_small,maxima_large, ...
               round(s1*Z1_sampling_fac), ...
               round(s2*Z2_sampling_fac));
           %
           nloc_small = sample_nlocs(points, maxima_small, round(Z1_sampling_fac*s1));
           nloc_large = sample_nlocs(points, maxima_large, round(Z2_sampling_fac*s2));           
       end
       %

disp(toc/60);
disp('clustering - done!')

h = figure();
ax = axes('parent',h);
set(ax, 'NextPlot', 'add');
x=out(:,1);
y=out(:,2);
ls=out(:,4);
is_bckg = 0==ls;
is_small = 1==ls;
is_large = 2==ls;
%
plot(ax, x(1==is_bckg), y(1==is_bckg), 'Marker', '.', 'MarkerSize', 5, 'LineStyle', 'none', 'color', rgb(127, 140, 141)); % gray
hold(ax,'on');
plot(ax, x(1==is_small), y(1==is_small), 'Marker', '.', 'MarkerSize', 5, 'LineStyle', 'none', 'color', rgb(53, 94, 59)); % green
hold(ax,'on');
plot(ax, x(1==is_large), y(1==is_large), 'Marker', '.', 'MarkerSize', 5, 'LineStyle', 'none', 'color', rgb(210, 4, 45)); % red
hold(ax,'on');
plot(ax, perim(:,1),perim(:,2), 'Marker', '.', 'MarkerSize', 2, 'LineStyle', 'none', 'color',rgb(44, 62, 80)); % dark grey
hold(ax,'off');
daspect(ax,[1 1 1]);
%
Name = 'gr_assisted_clustering_SL';
set(ax, 'box', 'on','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
grid(ax,'on');
saveas(h,[save_dir_name filesep Name '.fig']);
close(h);
%

%%%%%%%%%%%%%%%% ClusDoC-based - starts    
    class = out(:,3);
    datathr = out(:,1:2);

    display1 = true;
    printOutFig = true;
    maskVector = ones(size(class));

    DBSCANParams.epsilon = params.DBSCAN_epsilon;
    DBSCANParams.minPts = minPts1;
    DBSCANParams.SmoothingRad = params.DBSCAN_SmoothingRad;
    DBSCANParams.Cutoff = minPts1;
    DBSCANParams.DoStats = true;
    DBSCANParams.Outputfolder = save_dir_name;

        if display1 || ~printOutFig
            fig1 = figure();
            ax1 = axes('parent',fig1);
            set(ax1, 'NextPlot', 'add');
            plot(ax1, out(:,1), out(:,2), 'Marker', '.', 'MarkerSize', 5, 'LineStyle', 'none', 'color', rgb(127, 140, 141)); 
            axis image tight
        end

        SumofBigContour=[];
        SumofSmallContour=[];
        ClusterSmooth=cell(max(class),1);

        SumofContour = [];
        
        if isempty(class)
            % ?            
            SumofContour = {SumofBigContour, SumofSmallContour};                        
        else
            
            max_class = max(class);
            for i = 1:max_class

                xin = datathr(class == i,:); % Positions contained in the cluster i
                
                if contains_the_same_rows(xin(:,1:2))
                    ClusterSmooth{i,1} = [];
                    continue;
                end                    

                % assignin('base', 'xin', xin); % whatisthat

                if display1 || ~printOutFig
                    clusterColor = rand(1,3); % :)
                    plot(ax1,xin(:,1), xin(:,2),'Marker', '.', 'MarkerSize', 5,'LineStyle', 'none', 'color', clusterColor);
                end
                                                                                                                    
                [ClusImage,  Area, Circularity, Nb, contour, edges, Cutoff_point, Elongation] = Smoothing_fun4cluster(xin(:,1:2), DBSCANParams, false, false); % 0.1*max intensity 

                ClusterSmooth{i,1}.ClusterID = i;
                ClusterSmooth{i,1}.Points = xin(:,1:2);
                ClusterSmooth{i,1}.Area = Area;%
                ClusterSmooth{i,1}.Nb = Nb;%
                ClusterSmooth{i,1}.edges = edges;%
                ClusterSmooth{i,1}.Cutoff_point = Cutoff_point;
                ClusterSmooth{i,1}.Contour = contour;%
                ClusterSmooth{i,1}.Circularity = Circularity;%
                ClusterSmooth{i,1}.Elongation = Elongation;%
                ClusterSmooth{i,1}.Density_Nb_A = Nb/Area;%
                ClusterSmooth{i,1}.NInsideMask = sum(maskVector(class == i));
                ClusterSmooth{i,1}.NOutsideMask = sum(~maskVector(class == i));

                if Nb >= DBSCANParams.Cutoff
                    SumofBigContour = [SumofBigContour; contour; NaN NaN ];
                else
                    SumofSmallContour = [SumofSmallContour; contour; NaN NaN ];  
                end
                SumofContour={SumofBigContour, SumofSmallContour};

                % Plot the contour
                if display1 || ~printOutFig

                    if length(Nb) > DBSCANParams.Cutoff % Does this switch do anything?
                        plot(ax1, contour(:,1), contour(:,2), 'color', 'red');
                        set(ax1, 'box', 'on', 'XTickLabel', [], 'XTick', [], 'YTickLabel', [], 'YTick', []);
                    else
                        plot(ax1, contour(:,1), contour(:,2), 'color', rgb(44, 62, 80));
                    end

                end

            end
            
        end

        ClusterSmooth = ClusterSmooth(~cellfun('isempty', ClusterSmooth));

        disp(toc/60);
        disp('morph params - done!')

             Name = 'gr_assisted_clustering';
             set(ax1, 'box', 'on','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
             if printOutFig
                saveas(fig1,[DBSCANParams.Outputfolder filesep Name '.fig']);
                close(fig1);
             end

        Result = [];
        if DBSCANParams.DoStats

            Result.Number_Cluster = numel(ClusterSmooth);
            Result.Number(1) = mean(cell2mat(cellfun(@(x) x.Nb(x.Nb > DBSCANParams.Cutoff), ClusterSmooth, 'UniformOutput', false)));
            Result.Area(1) = mean(cell2mat(cellfun(@(x) x.Area(x.Nb > DBSCANParams.Cutoff), ClusterSmooth, 'UniformOutput', false)));
            Result.Mean_Circularity(1) = mean(cell2mat(cellfun(@(x) x.Circularity(x.Nb > DBSCANParams.Cutoff), ClusterSmooth, 'UniformOutput', false)));
            Result.Mean_Elongation(1) = mean(cell2mat(cellfun(@(x) x.Elongation(x.Nb > DBSCANParams.Cutoff), ClusterSmooth, 'UniformOutput', false)));
            %            
            Result.TotalNumber = size(datathr,1);
            Result.Percent_in_Cluster = sum(cell2mat(cellfun(@(x) x.Nb, ClusterSmooth, 'UniformOutput', false)))/size(datathr,1);
        end

        ClusterSmoothTable = ClusterSmooth; 
        save([DBSCANParams.Outputfolder filesep 'gr_Cluster_Result.mat'],'ClusterSmoothTable','Result','-v7.3');

        % save cluster info
        n1 = length(nloc_small);
        N1 = median(nloc_small);
        n2 = length(nloc_large);
        if 0~=n2
            N2 = median(nloc_large);
        else
            N2 = 0;
        end
        Ntot1_clusters = n1*N1; 
        Ntot2_clusters =  n2*N2; 
        Ntot_clusters = Ntot1_clusters + Ntot2_clusters;
        loc_percentage_in_clusters = Ntot_clusters/Ntot;
        rho1 = Ntot1_clusters/AREA;
        rho2 = Ntot2_clusters/AREA;
        rho1_clusters = n1/AREA;
        rho2_clusters = n2/AREA;
        caption = {'Z1','Z2','A1','A2','n1','n2','N1','N2','Ntot','Area','rho','Ntot1_clusters','Ntot2_clusters','Ntot_clusters','loc_percentage_in_clusters','rho1','rho2','rho1_clusters','rho2_clusters'};
        rec = [s1 s2 A1 A2 n1 n2 N1 N2 Ntot AREA rho Ntot1_clusters Ntot2_clusters Ntot_clusters loc_percentage_in_clusters rho1 rho2 rho1_clusters rho2_clusters];
        cell2csv([DBSCANParams.Outputfolder filesep 'gr_Clustering_2comp_stats.csv'],[caption; num2cell(rec)]);
        % save cluster info

        disp(toc/60);
        disp('+ ClusDoC params.. - done!')
%%%%%%%%%%%%%%%% ClusDoC-based - ends

end

% -------------FUNCTIONS------------------------------------

%-----------------------------------------------------------
function [out,v,ind,most_probable_v] = get_local_maxima(u,r)
%
z = u;
%
%smooth;
z = medfilt2(z,[3 3]);
z = imregionalmax(z);
% this needs compacting/cleaning  
z = imdilate(z,strel('disk',4));
zl = bwlabel(z);
s = regionprops(zl,'Centroid');
s = fliplr(round(vertcat(s.Centroid)));
out = s;
%
% remove multiples
d = squareform(pdist(s));
d = triu(d,1);
%
%close = d<(r/2) & d~=0; % ??
close = d<r & d~=0;

%
N = length(close);
Nmax = N^2;
multiples = zeros(Nmax,2);
cnt_multiples = 0;
for k=1:N
    for m=1:N
        if k > m, continue, end
        if close(k,m)
            cnt_multiples = cnt_multiples + 1;
            multiples(cnt_multiples,:) = [k m];
            %disp([k m d(k,m) r/2]);
        end
    end
end
multiples = multiples(1:cnt_multiples,:);

multiples_indices = unique(multiples);
good_indices  = setxor(1:N,multiples_indices);

multiple_maxima = s(multiples_indices,:);
good_maxima = s(good_indices,:);

fixed_maxima = [];
if ~isempty(multiple_maxima)
    idx = dbscan(multiple_maxima,r,2);
        for k=1:max(idx)
            coord = multiple_maxima(k==idx,:);
            fixed_maxima = [fixed_maxima; round([mean(coord(:,1)) mean(coord(:,2))])];
        end
end

out = cat(1,good_maxima,fixed_maxima);

v = u(sub2ind(size(u), out(:,1), out(:,2)) );

% remove too faint clusters (this is legacy & likely never has any effect)
[counts, edges] = histcounts(v(:));
max_bin_index = find(counts==max(counts)); 
most_probable_v = (edges(max_bin_index) + edges(max_bin_index + 1)) / 2;
most_probable_v = most_probable_v(1); % safety
good_indices = v > 0.25*most_probable_v; % tune?
%good_indices = v > 0.1*most_probable_v; 
%
out = out(good_indices,:);
v = v(good_indices,:);
ind = good_indices;

end

%----------------------------------------------------------- 
function z = imdraw(sx,sy,coord,v)
        z = zeros(sx,sy);           
        z(sub2ind(size(z), coord(:,1), coord(:,2))) = v;
end

%-----------------------------------------------------------
function mask = round_mask(r) % radius
    r  = round(r);
    if 0==rem(2,2), r = r+1; end
    mask = zeros(2*r+1);
    mask(r+1,r+1) = 1;
    mask = bwdist(mask)<=r;
end
%-----------------------------------------------------------
function n_locs = sample_nlocs(points, coord, r)    
   circle = round_mask(r);
   a = floor(length(circle)/2); % half-radius    
   pad = a + 4; % specify the width of padding
   %
   padded_points = padarray(points, [pad, pad], 0, 'both');
   % 
   x = coord(:,1) + pad;
   y = coord(:,2) + pad;
   n_locs = zeros(length(x),1);   
   for k=1:length(x)
      xc=x(k);
      yc=y(k);
      rx = (xc-a):(xc+a);
      ry = (yc-a):(yc+a);      
      %
      u = padded_points(rx,ry).*circle;
      n_locs(k) = sum(u(:));
   end    
end

%-----------------------------------------------------------
function  pruned = prune_small_size_maxima(sx,sy,maxima_small,maxima_large, large_r)

    circle = round_mask(large_r);
    a = floor(length(circle)/2); % half-radius
    %
    % prepare mask..
    excl = zeros(sx,sy);
    % pad it
    pad = a+4; % specify the width of padding
    padded_excl = padarray(excl, [pad, pad], 0, 'both');    

   for k=1:size(maxima_large,1)
      x = maxima_large(k,1) + pad;
      y = maxima_large(k,2) + pad;
          rx = (x-a):(x+a);
          ry = (y-a):(y+a);      
          padded_excl(rx,ry) = circle | padded_excl(rx,ry);  
   end   
   %
   excl_ind = [];
   for k=1:size(maxima_small,1)
      x = maxima_small(k,1) + pad;
      y = maxima_small(k,2) + pad;
      if padded_excl(x,y)
          excl_ind = [excl_ind; k]; % shouldn't be many
          padded_excl(x,y) = 0; % padded_excl(x,y) + 3; % :)
      end
   end
   %
   pruned = maxima_small;
   pruned(excl_ind,:) = [];

% debug
%    for k=1:size(pruned,1)
%       x = pruned(k,1) + pad;
%       y = pruned(k,2) + pad;
%       padded_excl(x,y) = 1;
%    end
%    icy_imshow(uint16(padded_excl),[num2str(size(maxima_small,1)) ' : ' num2str(size(pruned,1))]);

end

%-----------------------------------------------------------
function [out,perim]  = actual_point_clustering(points,maxima_small,maxima_large,small_r,large_r)
    %
    large_circle = round_mask(large_r);
    small_circle = round_mask(small_r);
    a_large = floor(length(large_circle)/2); % half-radius    
    a_small = floor(length(small_circle)/2); % half-radius 
    %
    pad = a_large+4; % specify the width of padding
    padded_points = padarray(points, [pad, pad], 0, 'both');
       
    % create a template
    N_small = size(maxima_small,1);
    N_large = size(maxima_large,1);
    % global label count
    label = 0;
    % scene is called "u"..
    u = zeros(size(padded_points));
    %
    %large    
    if ~isempty(maxima_large)
        for k=1:N_large
          x = maxima_large(k,1) + pad;
          y = maxima_large(k,2) + pad;
              rx = (x-a_large):(x+a_large);
              ry = (y-a_large):(y+a_large);      
              %
              label = label + 1;
              res = large_circle*label;
              u_k = u(rx,ry);
              res(u_k~=0) = u_k(u_k~=0);
              u(rx,ry) = res;
              %
              %label
        end
    end
    %small
    if ~isempty(maxima_small)    
        for k=1:N_small
          x = maxima_small(k,1) + pad;
          y = maxima_small(k,2) + pad;
              rx = (x-a_small):(x+a_small);
              ry = (y-a_small):(y+a_small);      
              %
              label = label + 1;
              res = small_circle*label;
              u_k = u(rx,ry);
              res(u_k~=0) = u_k(u_k~=0);
              u(rx,ry) = res;
              %
              %label
        end
    end
    %
    %icy_imshow(uint16(u));  
    %
    % in collaboration with ChatGPT
    [x,y] = find(padded_points == 1);    
    idx = u(sub2ind(size(u), x, y));
    %
    % convention on large-small - 0 bckg, 1 small 2 large
    ls = ones(size(idx));
    ls(0==idx) = 0; % bckg
    ls(0~=idx & idx<=N_large) = 2; %large
    %
    out = [(x-pad) (y-pad) idx ls];
    %
    [x_perim,y_perim] = find(1 == bwmorph(imgradient(u)>0,'thin',inf));
    perim = [x_perim-pad,y_perim-pad];
end

%-----------------------------------------------------------
function ret = contains_the_same_rows(x)
    ret = true;
    for k=2:size(x,1) 
        if 2 ~= sum(x(k,:)==x(1,:))
            ret = false;
        end
    end
end





























