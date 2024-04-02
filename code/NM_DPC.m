
% Please kindly cite the paper: Junyi Guan, Sheng Li, Xiongxiong He, and Jiajia Chen 
% "A novel clustering algorithm by adaptively merging sub- clusters based on the Normal-neighbor and Merging force"
% DOI:10.1007/s10044-021-00981-1
% Pattern Analysis and Applications,2021
% The code was rewritten by Junyi Guan in 2024

function [cl,NCLUST,time] = NM_DPC(data,k,epsilon_input,percent)
import Library.*
close all;
[n,~]=size(data);

nor_k_min = 5; %  minimum number of normal neighbors, default as 5
merge_setp1_thr = 0.8; % merge threshold of merging step 1, default as 0.8
merge_setp2_thr = 0.8; % merge threshold of merging step 2, default as 0.8
sharpness_gauss = 0.25; % the sharpness of a Gaussian distribution cluster
% Sharpness_gauss's theoretical value is approximately equal to 0.2349, but in our paper, we use 0.25 as default.
min_pts_thr = round(n/50); % if a cluster's points is less 'min_pts_thr', the points in the cluster are viewed as outliers to be assigned letter

tic;

y = pdist(data);
[~,N] =size(y);

dist = squareform(y);
rho = zeros(n,1);
[orderdists,orderindexs] = sort(dist,2);
kNN_index = orderindexs(:,1:k);
dc = mean(orderdists(:,round(n/100*percent)));

for i=1:n-1
    for j=i+1:n
        rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
        rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
    end
end
%% obtain normal-neigbor
bor_labels = zeros(n,5);
normal_k = zeros(n,1);

for i=1:n
    normal_k(i) = k;
    kNN_index_i = kNN_index(i,:);
    for j = (nor_k_min+2):k
        min_in_dist_cur_nei = min(dist(kNN_index_i(j),kNN_index_i(1:j-1)));
        %  min_in_dist: the minimum inner neighbor distance for current neighbor
        sum_in_dist = 0;
        for jj = 2:j-1
            min_in_dist_in_range = min(dist(kNN_index_i(jj),kNN_index_i(1:jj-1)));
            % min_in_dist_in_range: the minimum inner neighbor distances in the range of current neighbor
            sum_in_dist = min_in_dist_in_range + sum_in_dist;
        end
        aver_min_in_dist = sum_in_dist/(j-2);
        % aver_min_in_dist: the average minimum inner distance
        epsilon = min_in_dist_cur_nei / aver_min_in_dist;
        if epsilon > epsilon_input
            new_k = j-1;
            normal_k(i) = new_k;
            break
        end
    end
end
%% obtain sub-clusters
nneigh = zeros(n,1);
for i=1:n
    nor_k = normal_k(i);
    for cur_k=2:nor_k
        j = kNN_index(i,cur_k);
        if rho(i) < rho(j)
            nneigh (i) = j;
            break;
        end
    end
end
n_pk = 0;
for i = 1:n
    if nneigh(i) == 0
        n_pk = n_pk +1;        
    end 
end

peaks = find(nneigh==0);

sl=-1*ones(n,1); %% sl: sub-labels of points.
sl(peaks) = (1:n_pk); %% give unique sub-labels to density peaks.
[~,OrdRho]=sort(rho,'descend');

for i=1:n
    if (sl(OrdRho(i))==-1)
        sl(OrdRho(i))=sl(nneigh(OrdRho(i)));%% inherit sub-labels from NPN
    end
end

NCLUST = max(sl);

%% obtain border points
border_pts = [];
k_b = round(mean(normal_k)/2);

for i = 1:n
    nor_k = normal_k(i);
    if nor_k >= k_b
        bor_labels(i,1:k_b) = sl(kNN_index(i,1:k_b));
        diff_labels = unique(bor_labels(i,1:k_b));
    else
        bor_labels(i,1:nor_k) = sl(kNN_index(i,1:nor_k));
        bor_labels(i,nor_k+1:k_b) = sl(kNN_index(i,1));
        diff_labels = unique(bor_labels(i,1:k_b));
    end
    if length(diff_labels) > 1
        diff_labels = sort(diff_labels);
        border_pts = [[i diff_labels(1:2) rho(i) sl(i)]; border_pts];
    end
end

pk_rho_pts = [];

for i= 1:NCLUST
    subcluster_pts_num = length(find(sl == i));
    pk_rho_pts = [pk_rho_pts;[i rho(peaks(i)) subcluster_pts_num]];
end

border_pts = sortrows(border_pts, 2);
adj_sub_clusters = unique(border_pts(:,2:3),'rows');

%% merging step 1: find adjacent subclusters
[n_adj,~] = size(adj_sub_clusters);
adj_clusters_and_border_rho_step1 = [];

for i = 1:n_adj   
    border_pts_between_two_clusters = border_pts(border_pts(:,2) == adj_sub_clusters(i,1) & border_pts(:,3) == adj_sub_clusters(i,2),4:5);
    highest_border_rho_cl = border_pts_between_two_clusters(border_pts_between_two_clusters(:,1)== max(border_pts_between_two_clusters(:,1)),:);
    highest_border_rho = highest_border_rho_cl(1); 
    adj_clusters_and_border_rho_step1 = [adj_clusters_and_border_rho_step1; adj_sub_clusters(i,:) highest_border_rho];
end

%% merging step 1: obtain overlapping thickness between adjacent subclusters
overlapping_inf_step1 = [adj_clusters_and_border_rho_step1(:,1:2) zeros(n_adj,1)];

for i = 1:n_adj
    cl_1 = pk_rho_pts(overlapping_inf_step1(i,1),2:3);
    cl_2 = pk_rho_pts(overlapping_inf_step1(i,2),2:3);
    pts_num_of_cl1 = pk_rho_pts(overlapping_inf_step1(i,1),3);
    pts_num_of_cl2 = pk_rho_pts(overlapping_inf_step1(i,2),3); 
    if pts_num_of_cl1>pts_num_of_cl2
        big_rho = cl_1;
    else
        big_rho = cl_2;
    end  
    overlapping = adj_clusters_and_border_rho_step1(i,3)/big_rho(1);
    overlapping_inf_step1(i,3) = overlapping;
    
end


%% merging step 1: merge overlapping adjacent subclusters, according to parameter "merge_setp1_thr".
cl = sl;
for i = 1:n_adj
    cl_1 = overlapping_inf_step1(i,1);
    cl_2 = overlapping_inf_step1(i,2);
    if overlapping_inf_step1(i,3) > merge_setp1_thr
        cl(cl==cl_2) = cl_1;
        overlapping_inf_step1(overlapping_inf_step1(:,1) == cl_2,1) = cl_1;
        overlapping_inf_step1(overlapping_inf_step1(:,2) == cl_2,2) = cl_1;
    end
end

%% merging step 2: obtain the adjacent subclusters

adj_cls_step1= unique(overlapping_inf_step1(:,1:2),'rows');
adj_cls_step1 = adj_cls_step1(adj_cls_step1(:,1)~=adj_cls_step1(:,2),:);
adj_clusters_and_border_rho_step2 = [overlapping_inf_step1(:,1:2) adj_clusters_and_border_rho_step1(:,3)];
adj_clusters_and_border_highest_rho = [];
adj_clusters_and_border_mean_rho = [];

for i = 1:length(adj_cls_step1(:,1))
    adj_clusters = adj_cls_step1(i,:);
    idx = find(adj_clusters_and_border_rho_step2(:,1)== adj_clusters(1) & adj_clusters_and_border_rho_step2(:,2)== adj_clusters(2));
    boders1 = adj_clusters_and_border_rho_step2(idx,:);
    boders_rho_mean = mean(boders1(:,3));
    boders1 = boders1(boders1(:,3)==max(boders1(:,3)),:);
    adj_clusters_and_border_highest_rho = [adj_clusters_and_border_highest_rho;boders1];
    adj_clusters_and_border_mean_rho = [adj_clusters_and_border_mean_rho; boders1(:,1:2) boders_rho_mean];
end

labels = unique(cl);
NCLUST = length(labels);

for i= 1:NCLUST
    rhos_in_cl1 = rho(cl == labels(i));
    pts_nums_in_cl(i) = length(find(cl==labels(i)));
    max_rho_in_cl(i) = max(rhos_in_cl1);
    sharpness =  mean(abs(rhos_in_cl1 - mean(rhos_in_cl1)))/(max(rhos_in_cl1)-min(rho));
    merging_ability(i) = sharpness_gauss/sharpness;
end

%% merging step 2: obtain the overlapping thickness between adjacent subclusters

n_adj2 =length(adj_cls_step1(:,1));
overlapping_inf_step2 = [adj_cls_step1 zeros(n_adj2,1)];

for i = 1:length(adj_cls_step1(:,1))
    x1 = overlapping_inf_step2(i,1);
    x2 = overlapping_inf_step2(i,2);
    idx1 =  find(labels == x1);
    idx2 =  find(labels == x2);
    pts_num_of_cl1 = pts_nums_in_cl(idx1);
    pts_num_of_cl2 = pts_nums_in_cl(idx2);
    if pts_num_of_cl1>pts_num_of_cl2
        big_rho = max_rho_in_cl(idx1);
    else
        big_rho = max_rho_in_cl(idx2);
    end
    
    overlapping = adj_clusters_and_border_highest_rho(i,3)/big_rho;
    overlapping_inf_step2(i,3) = overlapping;
end

%% merging step 2: caculate the merging force and merging subclsuters to final clusters, according to parameter "merge_setp2_thr".
merging_force_inf = overlapping_inf_step2;

for i = 1:length(adj_cls_step1(:,1))
    label1 = overlapping_inf_step2(i,1);
    label2 = overlapping_inf_step2(i,2);
    overlapping = overlapping_inf_step2(i,3);
    
    merging_ability1 = merging_ability(labels == label1);
    merging_ability2 = merging_ability(labels == label2);
    
    avrg_merging_ability = (merging_ability1 + merging_ability2) / 2;
    merging_force_inf(i,3) = overlapping* avrg_merging_ability;
    
    if overlapping* avrg_merging_ability > merge_setp2_thr
        cl(cl== label2) = label1;
        overlapping_inf_step2(overlapping_inf_step2(:,1) == label2,1) = label1;
        overlapping_inf_step2(overlapping_inf_step2(:,2) == label2,2) = label1;
    end
end

%% updata information
labels = unique(cl);
NCLUST = length(labels);
for i = 1:NCLUST
    pts_nums_in_cl_new(i) = length(find(cl==labels(i)));
end

%% assign outlier points
outlier_cls_num = length(find(pts_nums_in_cl_new< min_pts_thr));
if outlier_cls_num>0 
    out_idx = find(pts_nums_in_cl_new< min_pts_thr);
    for i = 1:outlier_cls_num
        cl(cl == labels(out_idx(i))) = -1;
    end
    
    cl_new = cl;
    
    for i = 1:n
        if cl(i) == -1
            list = [(1:n)' dist(i,:)' cl];
            list = sortrows(list,2);
            for j = 1:n
                if list(j,end) > -1
                    cl_new(i,end) = list(j,end) ;
                    break;
                end
            end
            
        end
    end
    cl = cl_new;
    labels =  unique(cl);
    NCLUST = length(labels);
end

time = toc;

%% updata cluster labels
cl0 = cl;
for i = 1:length(labels)
    cl(cl0==labels(i)) = i;
end

end


