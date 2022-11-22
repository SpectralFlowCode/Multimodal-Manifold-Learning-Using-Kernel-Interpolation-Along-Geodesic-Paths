addpath(genpath('Utils'));

%% Load data and parameters
CM_GetParameters;
CM_GetData;

%% Figure 4
v2struct(SimParamsFigure2);
CM_AnalyzePair;

%% Figure 9 (in the suppplementary material)
v2struct(SimParamsFigure5);
CM_AnalyzePair;
