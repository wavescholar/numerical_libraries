function [G]=ImprovedFastGaussTransform(d,N,M,X,h,q,Y,p_max,K,ClusterIndex,ClusterCenter,ClusterRadii,r,epsilon)
%    Improved Fast Gauss Transform.
%
%     Computes and approximation $$\hat{G}(y_j)$$ to $$G(y_j)=\sum_{i=1}^{N} q_i
%     e^{\|x_i-y_j\|^2/h^2},\:\:j=1...M$$  such that
%     $$|\hat{G}(y_j)-G(y_j)| \leq Q \epsilon$$ , where
%     $$Q=\sum_{i=1}^{N}q_i$$.
%
%     C++ Implementation.
%
%     Loads ImprovedFastGaussTransform.dll
%
%    Implementation based on:
%
%     Fast computation of sums of Gaussians in high dimensions.  Vikas C. Raykar, C. Yang, R. Duraiswami, and N. Gumerov, CS-TR-4767, Department of computer science,University of Maryland, Collegepark.
%
% 
%% Input
%
%     * d                 --> data dimensionality.
%     * N                 --> number of source points.
%     * M                 --> number of target points.
%     * X                 --> d x N matrix of N source points in d dimensions.
%     * h                 --> the source scale or bandwidth.
%     * q                 --> 1 x N vector of the source strengths.
%     * Y                 --> d x M matrix of M target points in d dimensions.
%
%     * p_max        -->  maximum truncation number for the Taylor series.
%     * K                 -->  the  number of clusters.
%    * ClusterIndex --> 1 X N vector  the i th element is  the cluster   number  to which the i th point belongs. [ ClusterIndex[i] varies between 0 to K-1. ]
%    * ClusterCenter --> d x K matrix of K  cluster centers.
%    * ClusterRadii   --> 1 x K matrix of the radius of each cluster.
%     * r                   --> cutoff radius 
%     * epsilon         --> desired error
%
%% Ouput
%
%    * G                --> 1 x M vector of the Gauss Transform  evaluated at  each target point.
%
%% Signature
%
% Author: Vikas Chandrakant Raykar
% E-Mail: vikas@cs.umd.edu
% Date:  15 July 2005
%
%% See also
%
%  ImprovedFastGaussTransformChooseParameters,  ImprovedFastGaussTransformChooseTruncationNumber,  KCenterClustering,  example


