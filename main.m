% This code is for our paper: Gao, Chenqiang, et al. "Infrared small-dim target detection based on Markov random field guided noise modeling." Pattern Recognition 76 (2018): 463-475.
% If you use this code for your work, please cite this paper. 
% @article{Gao2018Infrared,
%    author = {Gao, Chenqiang and Wang, Lan and Xiao, Yongxing and Zhao, Qian and Meng, Deyu},
%    title = {Infrared small-dim target detection based on Markov random field guided noise modeling},
%    journal = {Pattern Recognition},
%    volume = {76},
%    number = {Supplement C},
%    pages = {463-475},
%    month = {2018/04/01/},
%    year = {2018}
% }

%% usage
% Some modules are written by C++, and you should recompile them if you run
% it with errors, Like as follows:
%>> mex update_Z.cpp
%>> mex reconstructImage.cpp

% readPath = './images'; % the path reading images
% savePath = './results'; % the path saving results, including background images and component images
% mog_IR_detection(readPath, savePath, temporal_step, patch, param, lr_prior, mog_prior)

close all;
clear all;
clc;
%% parameter setting
% patch parameter
temporal_step = 3;   % temporal sliding length and it is 3 frames at default.
patch.step = 5;
patch.size = 50; % 
patch.length = 3; % the number of frames for patching

%% model parameter
param.mog_k = 3; % the component number and it is 3 components at default
param.lr_init = 'SVD';
param.maxiter = 200;
param.initial_rank = 30;
param.tol = 1e-3;
lr_prior.a0 = 1e-6;
lr_prior.b0 = 1e-6;   

mog_prior.mu0 = 0;
mog_prior.c0 = 1e-6;
mog_prior.d0 = 1e-6;
mog_prior.alpha0 = 1e-6;
mog_prior.beta0 = 1e-6;

%% begin to process one image sequence
readPath = './images'; % the path reading images
savePath = './results'; % the path saving results, including background images and component images
mog_IR_detection(readPath, savePath, temporal_step, patch, param, lr_prior, mog_prior)