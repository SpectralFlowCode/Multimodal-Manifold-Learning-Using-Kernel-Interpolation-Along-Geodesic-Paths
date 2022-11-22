set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontsize',40);
rng(42);

%----output params----
if ~exist('OutputFigures')
    mkdir('OutputFigures')
end
OutputFolder='OutputFigures';
if ~exist('OutputTables')
    mkdir('OutputTables')
end
TablelsOutputFolder='OutputTables';
FigPreamble='NoisyMnist';


%----data params----
DataParams=struct;
DataParams.NumberOfDatapoints=1e3;
DataParams.MaxRotation=90;
DataParams.ShiftsRange=-1;
DataParams.RectSize=8;

%----simulation params----
SimParams.nReps=10;

%----kernels construnction params----
NormFac=0.2;
% 0<NormFac<1 - global scale factor that eqauls to NormFac*(median of the pairwise distances)
% NormFac>1 - local scale factor thar equals to the median of the distances of the NormFac nearest neighbours
UseMutualScaleFac='None';
% UseMutualScaleFac='Mean' - use the same scale factor for each kernel, the mutual scale factor is set to be their mean
% UseMutualScaleFac='Min' - use the same scale factor for each kernel, the mutual scale factor is set to be the lower scale factor
% UseMutualScaleFac='None' - use different scale factor for each kernel
KernelsParams=v2struct(NormFac,UseMutualScaleFac);

%----eigenvalues flow diagram params----
NumberOfEigenVals=10; %Number of eigenvalues calculated at each point on the geodesic
ntVec=200;% Number of points on the geodesic grid
TolFac=1e-3;% Tolerance factor for fixed rank approximation
Interpolator='Geodesic';
% Interpolator='Linear' - linear interpolation: (1-t)*K1+t*K2
% Interpolator='Geodesic' - geodesic interpolation: K1^(-1/2)*(K1^(1/2)*K2^(-1)*K1^(1/2))^t
% Interpolator='Harmonic' - Harmonic interpolation: K1^(1-t)*K2^t
EvfdParams=v2struct(NumberOfEigenVals,ntVec,TolFac);

%----alg params----
CMRParams=struct;
CMRParams.ntVec=100;
CMRParams.n=10; % the common spectrum
CMRParams.Ne=10; % at each point - number of components to check their commonality
CMRParams.acc_path=linspace(0.1,0.8,10);
CMRParams.Fast=true;
CMRParams.criterion='smoothness';




