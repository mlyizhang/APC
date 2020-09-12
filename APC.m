%% Experimental version of the proposed APC algorithm
%Wang Y, Wang D, Pang W, et al. A Systematic Density-based Clustering Method Using Anchor Points[J]. Neurocomputing, 2020.
%%
clc;clear;close all
addpath('D:\MEGAFile\work\evaluation', 'D:\MEGAFile\work\Complicate','D:\MEGAFile\work\UCI','D:\MEGAFile\work\drawGraph');
load ('compound.mat');% load your datasets
%% Parameters setting
%*************here gives a example of Compound***%%%%%%%%%
percent =1.9;%pct
mdelta=1.37;
Eps=0.85;
MinPts=3;
%%
shapeset =data;
distset = computeSimi(shapeset);% 
dc = computeDc(distset, percent);% 
rhos = getLocalDensity(distset, dc);%
[deltas, nneigh] = getDistanceToHigherDensity(distset, rhos);
min_delta=mdelta;
min_rho=0;
xstep=max(rhos)/10;%
showDG(rhos, deltas);% % decision graph
filter = (rhos > min_rho) & (deltas > min_delta);
cluster_num = sum(filter);%
ords = find(filter);%
cluster = zeros(size(rhos));% 
color = 1;
for i = 1:size(ords, 2)
    cluster(ords(i)) = color;% 
    color = color + 1;
end
%% 
[sorted_rhos, rords] = sort(rhos, 'descend');%
for i = 1:size(rords, 2) % 
    if cluster(rords(i)) == 0% 
        neigh_cluster = cluster(nneigh(rords(i)));
        assert(neigh_cluster ~= 0, 'neigh_cluster has not assign!');
        cluster(rords(i)) = neigh_cluster;% 
    end
end
% obtain micro-clusters
%% dividing density levels
xlim=get(gca,'Xlim');%
ylim=get(gca,'Ylim'); %
subspace={};% 
for i=1:(xlim(2)/xstep)
    subspace{i}=find(rhos([ords])>=(i-1)*xstep & rhos([ords])<=i*xstep);
end
ss=[];
for i=1:length(subspace)
    if isempty (subspace{i})
        ss(i)=0;
    else
        ss(i)=1;
    end
end
last1 = find(ss,1,'last');
first1=find(ss,1);
ss = ss(first1:last1); %
start=twozeroleast(ss);% 
if start==100
      fprintf('output results----idx');
   % 直接输出结果。
    idx=cluster';
Evaluation(label,idx);% output 
return;
end
idx=zeros(1,size(data,1));
silarm=[];%  
silar=[];
for i=first1:start+first1-1
temp=subspace{i};
silar=[silar temp];
end
for i=1:length(silar)
    a=find (cluster==silar(i));
    silarm=[silarm a];
end
idx([silarm])=200;%
%% outliers detection
DM=computeSimi(data([silarm],:)); %
%k=20
  k=1*size(DM,1)/3;
avgPneigh=avrgkneigbour(k,DM);%
x=1:length(avgPneigh);
newcl=[];
for i=1: length(avgPneigh)-1
    if avgPneigh(i+1)>=avgPneigh(i)*2
        if i>length(avgPneigh)/2 
            newcl=silarm([i+1:length(avgPneigh)]);
        else
            newcl=silarm([1:i]);% 
        end
    end
end
idx([newcl])=0;
datanum=[];%
for i=1:length(data(:,1))
    if idx(i)==0
       datanum(i)=i;
    end
end
hold off;
puri_datanum=find (datanum>0);% 
%% dbscan
distmat=computeSimi(data([puri_datanum],:));
Clust = DBSCAN(distmat,Eps,MinPts);
idx([puri_datanum])=Clust;
%  showdbResults(idx,data);
idx1=idx;
%% optimize dbscan results 
outlier_rep=[];
ordoutlier=[];% 
for i=1:length(ords)
     if idx(ords(i))==100
         ordoutlier=[ordoutlier 
             i];       
     end
end
ords([ordoutlier' silar]) =[]; %
for i=1:length(idx)
    if idx(i)==100 %
     [a,b]=  min( distset(i,[ords]));
     outlier_rep(i)=ords(b);
    end
end
for i=1:length(outlier_rep)
     if outlier_rep(i)>0
         idx(i)=idx(outlier_rep(i));
     end
end
ndata_cluster=[];%
numidx=unique (idx);
for ii=1:length(numidx)
     if numidx(ii)<100
         ndata_cluster(ii)=length(find(idx==ii));
     end
end
dbpeak=[];
dbrho=[];
dbnum=[];
points_idx=[];
dppeak_rho=[];
for jj=1:length(numidx)
    if numidx(jj)<=100
    points_idx=find(idx==jj);%
    [dbrho dbnum]=max (rhos(points_idx));
    dbpeak(jj)=points_idx(dbnum);%
    dppeak_rho(jj)=dbrho;% 
    end
end
dbpeak_simi=computeSimi(data(dbpeak,:));
[sdprho,sdprhonum]=sort(dppeak_rho, 'descend');
min_peak_cluster=min(find (sort(dppeak_rho, 'descend')<mean(dppeak_rho)));
if length(dppeak_rho)>2*(min_peak_cluster+0.5)
    for i=(min_peak_cluster+1):length(dppeak_rho)
    dbpeak_simi(i,i)=inf;
    bbb=sdprhonum(1:min_peak_cluster);
    [aa bb]=min (dbpeak_simi(sdprhonum(i),bbb));
    for jj=1:length(idx)
    if idx(jj)<=100
        idx(find(idx==i))=bbb(bb);
    end
end
end
end
%% report final results
fprintf('output results----idx\n');
idx=idx';
Evaluation(label,idx);