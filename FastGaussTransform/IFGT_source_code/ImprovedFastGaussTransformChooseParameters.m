function [K,p_max,r]=ImprovedFastGaussTransformChooseParameters(d,h,epsilon,Klimit)
%
%     Choose the parameters for the Improved Fast Gauss Transform.
%
%     C++ Implementation.
%
%     Loads ImprovedFastGaussTransformChooseParameters.dll
%
%
%    Implementation based on:
%   
%    Fast computation of sums of Gaussians in high dimensions. 
%    Vikas C. Raykar, C. Yang, R. Duraiswami, and N. Gumerov,
%    CS-TR-4767, Department of computer science,
%    University of Maryland, Collegepark.
%
%% Input
%
% * d...dimension of the points.
% * h... the source bandwidth.
% * epsilon...the desired error.   
% * Klimit ...upper limit on the number of clusters, Klimit. 
%
%  Note :  [ Use roughly Klimit=round(40*sqrt(d)/h) ] 
%
%% Ouput
%
% * K...number of clusters.
% * p_max...maximum truncation number for the Taylor series.
% * r...source cutoff radius.
%
%% Signature
%
% Author: Vikas Chandrakant Raykar
% E-Mail: vikas@cs.umd.edu
% Date: August 22, 2005
%

