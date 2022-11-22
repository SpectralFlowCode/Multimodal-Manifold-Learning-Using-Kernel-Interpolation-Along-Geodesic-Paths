% clc;
% clear all;
% close all;
addpath(genpath('Utils'));

%% Load data and parameters
ENose_GetParameters;
ENose_GetData;

%% Figure 5
v2struct(SimParamsFigure3);
ENose_AnalyzePair;


