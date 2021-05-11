function gather_SMLM_quantifiers(path_to_data,Ch1_name,Ch2_name,Ch1_label,Ch2_label)

Ch1_ClusDoC_dir = [path_to_data filesep Ch1_name '_channel_1_ClusDoC_Results'];
Ch2_ClusDoC_dir = [path_to_data filesep Ch2_name '_channel_2_ClusDoC_Results'];
Ch1_MIiSR_dir = [path_to_data filesep Ch1_name '_channel_1_MIiSR_Results'];
Ch2_MIiSR_dir = [path_to_data filesep Ch2_name '_channel_2_MIiSR_Results'];

% ROIs and CLUSTERS
% DoC
DBSCAN_DoC_data_fname = [Ch1_ClusDoC_dir filesep 'DoC' filesep 'DBSCAN Results' filesep 'DBSCAN Clus-DoC Results.mat'];
%DBSCAN_DoC_data_fname = [Ch1_ClusDoC_dir filesep 'DoC' filesep 'DBSCAN_Results' filesep 'DBSCAN_Clus-DoC_Results.mat'];

% ROIs
% ClusDoC Ripley
ClusDoC_Ripley_dirname_Ch1 = [Ch1_ClusDoC_dir filesep 'RipleyK' filesep 'RipleyK_results'];
ClusDoC_Ripley_dirname_Ch2 = [Ch2_ClusDoC_dir filesep 'RipleyK' filesep 'RipleyK_results'];

% ROIs
% SAA
MIiSR_SAAdata_fname = [Ch1_MIiSR_dir filesep 'SAA' filesep 'SAAstruct.mat'];


% ROIs
% Spatial data (cross RDF, cross Ripley)
MIiSR_SpatialStats_Ch1_data_fname = [Ch1_MIiSR_dir filesep 'SpatialStats' filesep 'SpatialData.mat'];
MIiSR_SpatialStats_Ch2_data_fname = [Ch2_MIiSR_dir filesep 'SpatialStats' filesep 'SpatialData.mat'];

tic

% ClusDoC Ripley
%
%Ch1
fnames = dir([ClusDoC_Ripley_dirname_Ch1 filesep '*.xls']);
fnames = {fnames.name};
RipelyKx = [];
RipelyK_Ch1_data = [];
roi_index_RipelyK_Ch1 = [];
for k=1:numel(fnames)
    if contains(fnames{k},'_ROI')
        [~,raw] = xlsread([ClusDoC_Ripley_dirname_Ch1 filesep fnames{k}]);
        x = raw(:,1);
        c = raw(:,2);
        x = x(2:size(raw,1));
        c = c(2:size(raw,1));
        if isempty(RipelyKx)
            RipelyKx = str2num(char(x));
        end
        c = str2num(char(c));        
        RipelyK_Ch1_data = [RipelyK_Ch1_data c];        
        %
        s = strsplit(fnames{k},'_ROI_');
        s = s{2};
        s=strrep(s,'.xls','');
        ind_k = str2num(s);
        roi_index_RipelyK_Ch1 = [roi_index_RipelyK_Ch1; ind_k];
    end
end
%
%Ch2
fnames = dir([ClusDoC_Ripley_dirname_Ch2 filesep '*.xls']);
fnames = {fnames.name};
RipelyKx = [];
RipelyK_Ch2_data = [];
roi_index_RipelyK_Ch2 = [];
for k=1:numel(fnames)
    if contains(fnames{k},'_ROI')
        [~,raw] = xlsread([ClusDoC_Ripley_dirname_Ch2 filesep fnames{k}]);
        x = raw(:,1);
        c = raw(:,2);
        x = x(2:size(raw,1));
        c = c(2:size(raw,1));
        if isempty(RipelyKx)
            RipelyKx = str2num(char(x));
        end
        c = str2num(char(c));
        RipelyK_Ch2_data = [RipelyK_Ch2_data c];
        %
        s = strsplit(fnames{k},'_ROI_');
        s = s{2};
        s=strrep(s,'.xls','');
        ind_k = str2num(s);
        roi_index_RipelyK_Ch2 = [roi_index_RipelyK_Ch2; ind_k];        
    end
end
% x1 = mean(RipelyK_Ch1_data,2);
% x2 = mean(RipelyK_Ch2_data,2);
% plot(RipelyKx,x1,'k.-',RipelyKx,x2,'r.-');

data.RipelyKx = RipelyKx;
data.RipelyK_Ch1_data = RipelyK_Ch1_data;
data.roi_index_RipelyK_Ch1 = roi_index_RipelyK_Ch1;
data.RipelyK_Ch2_data = RipelyK_Ch2_data;
data.roi_index_RipelyK_Ch2 = roi_index_RipelyK_Ch2;

% DBSCAN 
load(DBSCAN_DoC_data_fname);
%
DoC_DBSCAN_clusters_ch1 = [];
for k=1:numel(ClusterSmoothTableCh1)
    clusters_k = ClusterSmoothTableCh1{k};
    for c=1:numel(clusters_k)
        Area = clusters_k{c}.Area;
        Circularity = clusters_k{c}.Circularity;
        MeanDoC = clusters_k{c}.MeanDoC;
        N = clusters_k{c}.NInsideMask;
        DoC_DBSCAN_clusters_ch1 = [DoC_DBSCAN_clusters_ch1; [Area Circularity MeanDoC N]];
    end
end

DoC_DBSCAN_clusters_ch2 = [];
for k=1:numel(ClusterSmoothTableCh2)
    clusters_k = ClusterSmoothTableCh2{k};
    for c=1:numel(clusters_k)
        Area = clusters_k{c}.Area;
        Circularity = clusters_k{c}.Circularity;
        MeanDoC = clusters_k{c}.MeanDoC;
        N = clusters_k{c}.NInsideMask;        
        DoC_DBSCAN_clusters_ch2 = [DoC_DBSCAN_clusters_ch2; [Area Circularity MeanDoC N]];
    end
end
%
data.DoC_DBSCAN_clusters_ch1 = DoC_DBSCAN_clusters_ch1;
data.DoC_DBSCAN_clusters_ch2 = DoC_DBSCAN_clusters_ch2;

% ROIs
% SAA
load(MIiSR_SAAdata_fname);
data.SAAstruct = SAAstruct;

% ROIs
% Spatial data (cross RDF, cross Ripley)
load(MIiSR_SpatialStats_Ch1_data_fname);
data.SpatialData_Ch1 = SpatialData;
load(MIiSR_SpatialStats_Ch2_data_fname);
data.SpatialData_Ch2 = SpatialData;

data.Ch1_label = Ch1_label;
data.Ch2_label = Ch2_label;

save([path_to_data filesep Ch1_name],'data');
toc/60

end
