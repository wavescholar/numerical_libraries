% Script the demonstrate the use of IFGT

clear all;
close all;
clear functions;
clc;

disp('---------------------------------------------');
disp(sprintf('Example to demonstrate the use of IFGT'));
disp('---------------------------------------------');

% the data dimensionality

d=2;

disp(sprintf('Dimensionality d=%d\n',d));

% the number of sources

N=5000;

disp(sprintf('Number of source points N=%d\n',N));

% the number of targets

M=5000;

disp(sprintf('Number of target points M=%d\n',N));

% the source points
% d x N matrix of N source points in d dimensions.
% Scale the data

G=30; 
m=rand(d,G);
v=0.02*ones(1,G);
[X]=generate_multiple_gaussians(N,G,m,v,d);

for j=1:d
    shift=min(X(j,:));
    X_shifted(j,:)=X(j,:)-shift;
    scale=1/max(X_shifted(j,:));
    X_shifted_scaled(j,:)=X_shifted(j,:)*scale;
end

X= X_shifted_scaled;

% the target points
% d x M matrix of M source points in d dimensions.

Y=rand(d,M);

% the source weights
% 1 x N row vector

q=rand(1,N);  

% the bandwidth

h=sqrt(2)*0.2*sqrt(d);

disp(sprintf('Bandwidth h=%f\n',h));

% the desired error

epsil=1e-3;     

disp(sprintf('Target error epsilon=%e\n',epsil));

% The upperlimit on the number of clusters.
% If the number of clusters chosen is equal to Klimit then increase this
% value

Klimit=round(0.2*100/h);

disp(sprintf('Klimit=%d\n',Klimit));

disp(sprintf('Press any key to continue...\n'));

pause

disp('---------------------------------------------');
disp(sprintf('Choosing the IFGT parameters\n'));
disp('---------------------------------------------');

% Choose the parameters
%
% K          -- number of clusters
% p_max -- maximum truncation number
% r          -- cutoff radius

to=clock;
[K,p_max,r]=ImprovedFastGaussTransformChooseParameters(d,h,epsil,Klimit);
parameters_time=etime(clock,to);

disp(sprintf('Number of clusters K=%d\n',double(K)));
disp(sprintf('Maximum truncation number p_max=%d\n',double(p_max)));
disp(sprintf('Cutoff radius r=%f\n',r));
disp(sprintf('Time taken=%f secs\n',parameters_time));

disp(sprintf('Press any key to continue...\n'));

pause

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
title('Results of the k-center clustering procedure');

disp(sprintf('Press any key to continue...\n'));

pause

disp('---------------------------------------------');
disp(sprintf('Updating the truncation number\n'));
disp('---------------------------------------------');

to=clock;
[p_max]=ImprovedFastGaussTransformChooseTruncationNumber(d,h,epsil,rx);
trunc_time=etime(clock,to);

disp(sprintf('Updated Maximum Truncation Number=%d\n',double(p_max)));
disp(sprintf('Time taken=%f secs\n',trunc_time));


disp(sprintf('Press any key to continue...\n'));

pause

disp('---------------------------------------------');
disp(sprintf('Running the  IFGT\n'));
disp('---------------------------------------------');

to=clock;
[G_IFGT]=ImprovedFastGaussTransform(d,N,M,X,h,q,Y,double(p_max),double(K),ClusterIndex,ClusterCenter,ClusterRadii,r,epsil);
IFGT_time=etime(clock,to);

disp(sprintf('Time taken=%f secs\n',IFGT_time));

disp(sprintf('Press any key to continue...\n'));


pause

disp('---------------------------------------------');
disp(sprintf('Running the direct method.\n'));
disp('---------------------------------------------');

to=clock;
[G_direct]=GaussTransform(d,N,M,X,h,q,Y);
GT_time=etime(clock,to);

disp(sprintf('Time taken=%f secs\n',GT_time));

disp('---------------------------------------------');
disp(sprintf('Summary\n'));
disp('---------------------------------------------');

IFGT_total_time=parameters_time+clustering_time+trunc_time+IFGT_time;
IFGT_err=max(abs((G_direct-G_IFGT)))/sum(q);

disp(sprintf('Direct computation takes %f secs\n',GT_time));
disp(sprintf('IFGT takes %f secs Speedup=%f\n',IFGT_total_time,GT_time/IFGT_total_time));
disp(sprintf('Actual error for IFGT is %e Target was %e \n',IFGT_err,epsil));
disp('---------------------------------------------');



