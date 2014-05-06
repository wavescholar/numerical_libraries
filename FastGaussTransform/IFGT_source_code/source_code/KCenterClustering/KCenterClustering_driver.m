% Script the demonstrate the use of KCenterClustering

clear all;
close all;
clear functions;
clc;

disp('---------------------------------------------');
disp(sprintf('Example to demonstrate the use of  KcenterClustering'));
disp('---------------------------------------------');

% the data dimensionality

d=2;

disp(sprintf('Dimensionality d=%d\n',d));

% the number of sources

N=5000;

disp(sprintf('Number of source points N=%d\n',N));

% the source points
% d x N matrix of N source points in d dimensions.

X=rand(d,N);

% Number of clusters

K=20;

disp(sprintf('Number of clusters K=%d\n',K));

disp('---------------------------------------------');
disp(sprintf('Running the k-center clustering\n'));
disp('---------------------------------------------');

% k-center clustering

to=clock;
[rx,ClusterIndex,ClusterCenter,NumPoints,ClusterRadii]=KCenterClustering(d,N,X,double(K));
clustering_time=etime(clock,to);

disp(sprintf('Maximum cluster radius=%f\n',rx));
disp(sprintf('Time taken=%f secs\n',clustering_time));

%Pretty plot of the results of the clustering procedure. 
%Plots only for two and three dimensions.

plot_clusters(N,d,X,K,ClusterIndex,ClusterCenter);


