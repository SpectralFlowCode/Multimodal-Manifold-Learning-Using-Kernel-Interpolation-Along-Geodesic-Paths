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

%----data params----
DatasetsParams=struct;
DatasetsParams.DownloadLink='https://archive.ics.uci.edu/ml/machine-learning-databases/00447/data.zip';

%----alg params----
CMRParams=struct;
CMRParams.ntVec=100;
CMRParams.n=20; % the common spectrum
CMRParams.Ne=5; %at each point - number of components to check their commonality
CMRParams.acc_path=linspace(0.2,0.6,10);
CMRParams.Fast=true;
CMRParams.criterion='smoothness';

%----smoothness params----
SmoothnessParams=struct;
SmoothnessParams.PreProcess='none';
SmoothnessParams.NormFac=0.5;
SmoothnessParams.l=10;
SmoothnessParams.UseGS=true;

%----params for fig4----
[NormFac1,NormFac2,NumberOfEigenVals,ntVec]=deal(1,1,10,200);
[TargetInd,SensorInd1,SensorInd2]=deal(2,2,12);
PreProcess='none';
noise1='sine';As=0.1;NoiseParams1=struct('As',As,'f',200);noise2=noise1;NoiseParams2=struct('As',As,'f',400);
th=1e-3;
FigPreamble="CM_Figure4";
figCMRParams=CMRParams;
figCMRParams.Ne=2;
figCMRParams.n=10;
SimParamsFigure2=v2struct(FigPreamble,NormFac1,NormFac2,NumberOfEigenVals,ntVec,TargetInd,SensorInd1,SensorInd2,...
                    PreProcess,noise1,noise2,NoiseParams1,NoiseParams2,th,figCMRParams);

%----params for fig9----
[NormFac1,NormFac2,NumberOfEigenVals,ntVec]=deal(0.5,0.5,10,200);
[TargetInd,SensorInd1,SensorInd2]=deal(1,5,11);
PreProcess='median';
noise1='none';NoiseParams1=struct;noise2=noise1;NoiseParams2=struct;
th=5e-3;
FigPreamble="CM_Figure9";
figCMRParams=CMRParams;
figCMRParams.Ne=5;
SimParamsFigure5=v2struct(FigPreamble,NormFac1,NormFac2,NumberOfEigenVals,ntVec,TargetInd,SensorInd1,SensorInd2,...
                    PreProcess,noise1,noise2,NoiseParams1,NoiseParams2,th,figCMRParams);

                
