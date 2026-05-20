% Run this first time to add paths
clear
clear all classes
close all

% Add path but include the full directory name
addpath(fullfile(pwd,'./classes'));
% addpath(fullfile(pwd,'./analytic-cont'))
addpath(fullfile(pwd,'./functions'));
% Add the data folder
addpath(fullfile(pwd,'./../finite-stokeswave-data/'));