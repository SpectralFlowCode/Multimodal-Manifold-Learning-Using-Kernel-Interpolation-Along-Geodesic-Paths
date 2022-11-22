set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontsize',40);
rng(1)

%----output params----
if ~exist('OutputFigures')
    mkdir('OutputFigures')
end
OutputFolder='OutputFigures';
if ~exist('OutputTables')
    mkdir('OutputTables')
end
TablelsOutputFolder='OutputTables';

%----smoothness params----
SmoothnessParams=struct;
SmoothnessParams.PreProcess='zscore';
SmoothnessParams.NormFac=0.5;
SmoothnessParams.l=5;
SmoothnessParams.UseGS=true;

%----data params----
DatasetsParams=struct;
DatasetsParams.DownloadLink='http://archive.ics.uci.edu/ml/machine-learning-databases/00270/driftdataset.zip';
DatasetsParams.SensorsListPerGas={[2,7,11,14],[2,7,11,14],[6,11,14],[2,7,11,14],[2,7,11,14],[2,7,11,15]};

%----alg params----
CMRParams=struct;
CMRParams.ntVec=100;
CMRParams.n=20; % the common spectrum
CMRParams.Ne=5; %at each point - number of components to check their commonality
CMRParams.acc_path=linspace(0,1,10);
CMRParams.acc_path=0.5;
CMRParams.Fast=true;
CMRParams.criterion='smoothness';

%----params for fig4----
[NormFac1,NormFac2,NumberOfEigenVals,ntVec]=deal(1,1,10,200);
[TargetInd,SensorInd1,SensorInd2]=deal(4,5,11);
common_insets={[0,3,4],[0.4,2,4],[1,4,7]};
noise_insets={[0,2,5],[1,2,3]};
PreProcess='zscore';
PolyNormFac=0.5;
FigPreamble="ENose";
SimParamsFigure3=v2struct(FigPreamble,NormFac1,NormFac2,NumberOfEigenVals,ntVec,TargetInd,SensorInd1,SensorInd2,...
                    common_insets,noise_insets,PreProcess,PolyNormFac);
                

