addpath(genpath('Utils'));

%% 2D Tori
% Commmon is the polodial angle (Figure 2 and Figure 3 in the paper)
CommonPolodial_GetParameters;
Simulations_GetData;
Simulations_AnalyzePair;

% Commmon is the torodial angle (Figure 24 and Figure 25 in the SM)
CommonTorodial_GetParameters;
Simulations_GetData;
Simulations_AnalyzePair;

%% 2D flat manifold 
% Moderate SNR (Figure 14 and Figure 15 in the SM)
TwoDFlat_GetParameters;
Simulations_GetData;
Simulations_AnalyzePair;

    