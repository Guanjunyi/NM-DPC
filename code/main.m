% Please kindly cite the paper: Junyi Guan, Sheng Li, Xiongxiong He, and Jiajia Chen 
% "A novel clustering algorithm by adaptively merging sub- clusters based on the Normal-neighbor and Merging force"
% DOI:10.1007/s10044-021-00981-1
% Pattern Analysis and Applications,2021
% The code was rewritten by Junyi Guan in 2024

clear;close all;clc;
%% load dataset
load data/Flame
data_with_lable = Flame;
%% deduplicate data
data_x = unique(data_with_lable,'rows');
if size(data_x,1) ~= size(data_with_lable,1)
    data_with_lable = data_x;
end
lable = data_with_lable(:,end);
data = data_with_lable(:,1:end-1);
data=(data-min(data))./(max(data)-min(data));
data(isnan(data))=0;
%% NM-DPC clustering
k = 20;
epsilon_input = 2;
percent = 2;
[CL,NC,runtime] =  NM_DPC(data,k,epsilon_input,percent);
%% evaluation
[AMI,ARI,FMI] = Evaluation(CL,lable);
%% show result
resultshow(data,CL);
%% clustering result
result = struct;
result.NC = NC;
result.AMI = AMI;
result.ARI = ARI;
result.FMI = FMI;
result.runtime = runtime;
result