function [p_max]=ImprovedFastGaussTransformChooseTruncationNumber(d,h,epsilon,rx)
%
%     Choose the parameters for the Improved Fast Gauss Transform.
%
%     Given the cluster radius returns the truncation number required.
%
%     C++ Implementation.
%
%     Loads ImprovedFastGaussTransformChooseTruncationNumber.dll
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
% * rx...maximum cluster radius
%
%% Ouput
%
% * p_max...maximum truncation number for the Taylor series.
%
%% Signature
%
% Author: Vikas Chandrakant Raykar
% E-Mail: vikas@cs.umd.edu
% Date: August 22, 2005
%
