% Script the demonstrate the use of GaussTransform

clear all;
close all;
clear functions;
clc;

disp('---------------------------------------------');
disp(sprintf('Example to demonstrate the use of Gauss Transform'));
disp('---------------------------------------------');

% the data dimensionality

d=2;

disp(sprintf('Dimensionality d=%d\n',d));

% the number of sources

N=100;

disp(sprintf('Number of source points N=%d\n',N));

% the number of targets

M=100;

disp(sprintf('Number of target points M=%d\n',N));

% the source points
% d x N matrix of N source points in d dimensions.

X=rand(d,N);

% the target points
% d x M matrix of M source points in d dimensions.

Y=rand(d,M);

% the source weights
% 1 x N row vector

q=rand(1,N);  

% the bandwidth

h=sqrt(2)*0.2*sqrt(d);

disp(sprintf('Bandwidth h=%f\n',h));

disp('---------------------------------------------');
disp(sprintf('Running the Gauss Transform.\n'));
disp('---------------------------------------------');

to=clock;
[G_direct]=GaussTransform(d,N,M,X,h,q,Y);
GT_time=etime(clock,to);

disp(sprintf('Time taken=%f secs\n',GT_time));


