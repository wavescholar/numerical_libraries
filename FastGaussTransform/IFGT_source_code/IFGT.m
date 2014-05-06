function [G]=IFGT(d,N,M,X,h,q,Y,epsil)
%
%     Fast computation of the Gauss Transform.
%
%     Computes and approximation $$\hat{G}(y_j)$$ to $$G(y_j)=\sum_{i=1}^{N} q_i
%     e^{\|x_i-y_j\|^2/h^2},\:\:j=1...M$$  such that
%     $$|\hat{G}(y_j)-G(y_j)| \leq Q \epsilon$$ , where
%     $$Q=\sum_{i=1}^{N}q_i$$.
%
%     C++ Implementation.
%
%
%% Input
%
% * d                 --> data dimensionality.
% * N                 --> number of source points.
% * M                 --> number of target points.
% * X                 --> d x N matrix of N source points in d dimensions.
% * h                 --> source bandwidth or scale.
% * q                 --> 1 x N vector of the source strengths.
% * Y                 --> d x M matrix of M target points in d dimensions.
% * epsil         --> desired error
%
%% Ouput
%
% * G                --> 1 x M vector of the Gauss Transform evaluated at  each target point.
%
%% Signature
%
% Author: Vikas Chandrakant Raykar
% E-Mail: vikas@cs.umd.edu
% Date:  08 July 2005
%
%% See also
%
%  ImprovedFastGaussTransformChooseParameters,  ImprovedFastGaussTransformChooseTruncationNumber,  KCenterClustering,  ImprovedFastGaussTransform
%


Klimit=round(0.2*sqrt(d)*100/h);
% The upperlimit on the number of clusters.
% If the number of clusters chosen is equal to Klimit then increase this
% value


disp(sprintf('Klimit=%d\n',Klimit));

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

disp('---------------------------------------------');
disp(sprintf('Running the k-center clustering\n'));
disp('---------------------------------------------');

% k-center clustering

to=clock;
[rx,ClusterIndex,ClusterCenter,NumPoints,ClusterRadii]=KCenterClustering(d,N,X,double(K));
clustering_time=etime(clock,to);

disp(sprintf('Maximum cluster radius=%f\n',rx));
disp(sprintf('Time taken=%f secs\n',clustering_time));

disp('---------------------------------------------');
disp(sprintf('Updating the truncation number\n'));
disp('---------------------------------------------');

to=clock;
[p_max]=ImprovedFastGaussTransformChooseTruncationNumber(d,h,epsil,rx);
trunc_time=etime(clock,to);

disp(sprintf('Updated Maximum Truncation Number=%d\n',double(p_max)));
disp(sprintf('Time taken=%f secs\n',trunc_time));

disp('---------------------------------------------');
disp(sprintf('Running the  IFGT\n'));
disp('---------------------------------------------');

to=clock;
[G]=ImprovedFastGaussTransform(d,N,M,X,h,q,Y,double(p_max),double(K),ClusterIndex,ClusterCenter,ClusterRadii,r,epsil);
IFGT_time=etime(clock,to);

disp(sprintf('Time taken=%f secs\n',IFGT_time));

