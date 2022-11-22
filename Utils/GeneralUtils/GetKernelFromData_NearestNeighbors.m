function [ Kx] = ...
    GetKernelFromData_NearestNeighbors(X, d, params)
% NCCA: Nonparametric canonical correlation analysis

%%%%%
% Unlike NCCA, in NCCA2 the decomposition is done accurately (and less efficiently ) via SVD decompoistion (see line 180)
%%%%%


%
% Input:
% X,Y - paired training examples (rows) of features (columns).
% d - dimension of the output transformed features.
% X_unpaired, Y_unpaired [optional] - additional unpaired examples. The
%       numbers of unpaired X and Y examples does not need to be identical.
% params [optional] - structure with algorithm paramters:
%       - hx,hy - bandwidth parameters for the KDEs of views 1,2. Default
%       is 0.5.
%       - nnx,nny - number of nearest neighbors for the KDEs of views 1,2.
%       Default is 20.
%       - randSVDiters - number of iterations for random SVD algorithm.
%       Set higher for better accuracy, default is 20.
%       - randSVDblock - block size for random SVD algorithm. Set higher
%       for better accuracy, default is d+10;
%       - doublyStIters - number of iterations for doubly stichastic
%       normalization. Set higher for better accuracy, default is 15.
%
%
% Output:
% Kx, Ky - sparse kernels
%
% This code is based on the paper:
%       Tomer Michaeli, Weiran Wang and Karen Livescu,
%       "Nonparametric Canonical Correlation Analysis",
%       International Conference on Machine Learning (ICML 2016)
%
% Written by Tomer Michaeli (Technion - IIT) and Weiran Wang (TTIC).
%
% This Matlab code is distributed only for academic research purposes.
% For other purposes, please contact Tomer Michaeli
% mail: tomer.m@ee.technion.ac.il
%
% Version 2.0, January 2017.


%% Check inputs and set default values to the unspecified parameters
narginchk(3,6);

if ~exist('params','var')
    params = [];
end
if ~isfield(params,'hx')
    hx = 0.5;
else
    hx = params.hx;
end
if ~isfield(params,'nnx')
    nnx = 20;
else
    nnx = params.nnx;
end


N = size(X,1); % Number of training points


%% Normalize data
meanX = mean(X); 
X = bsxfun(@minus, X, meanX);

meanSqX = sqrt(mean(sum(X.^2,2))); 
X = bsxfun(@times, X, 1./meanSqX);


%% Compute NNs
% fprintf('Computing nearest neighbors ...'); tic;

[idxs_X,dists_X] = knnsearch(X(1:N,:), X, 'K', nnx);
dists_X = dists_X.^2'; idxs_X = idxs_X';
clear X


%% Compute weight matrices for the two views
colInd = ones(nnx, 1) * (1:(N));
Dx = sparse(double(idxs_X(:)), colInd(:), exp(-0.5*dists_X(:)/hx^2), N, N)';
clear colInd idxs_X idxs_Y dists_X dists_Y

Kx=0.5*(full(Dx)+full(Dx)');

end