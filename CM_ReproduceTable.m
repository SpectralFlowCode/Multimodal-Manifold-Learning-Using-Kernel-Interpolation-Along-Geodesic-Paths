% clc;
% clear all;
% close all;
addpath(genpath('Utils'));

%% Load data and parameters
CM_GetParameters;
CM_GetData;

%% Evaluate all pairs
disp(sprintf('Generating condition monitoring objective results'))
TargetInd=2;SensorsList=[1,2,3,11,12,15];

Target=ConditionNames{TargetInd};
IndsToProcess=RelevantInds(460:end);
GT=Conditions(Target);GT=GT(IndsToProcess);
cX=GT-mean(GT);

Lmax=SmoothnessParams.l;
ScoresPerNoise=[];
reverseStr = '';
for noise=["sine","saw"]
    msg = sprintf('Evaluating %s wave', noise);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    switch noise
        case 'sine'
            As=0.25;NoiseParams1=struct('As',As,'f',200);NoiseParams2=struct('As',As,'f',400);
        case 'saw'
            As=0.25;NoiseParams1=struct('Apd',As,'fd',5);NoiseParams2=struct('Apd',As,'fd',8);
    end
    
    KernelsList1={};DatapointsList1={};
    for SensorInd=SensorsList
        SensorToUse=SensorsToLoad{SensorInd};
        [ColumnStochasticK,K,NoisyData] =GetKernelFromNoisySensor(Samples,SensorToUse,SmoothnessParams.PreProcess,IndsToProcess,SmoothnessParams.NormFac,noise,NoiseParams1,false);
        KernelsList1=[KernelsList1,K];
        DatapointsList1=[DatapointsList1,NoisyData];
    end
    KernelsList2={};DatapointsList2={};
    for SensorInd=SensorsList
        SensorToUse=SensorsToLoad{SensorInd};
        [ColumnStochasticK,K,NoisyData] =GetKernelFromNoisySensor(Samples,SensorToUse,SmoothnessParams.PreProcess,IndsToProcess,SmoothnessParams.NormFac,noise,NoiseParams2,false);
        KernelsList2=[KernelsList2,K];
        DatapointsList2=[DatapointsList2,NoisyData];
    end
    
    SmoothnessLinear=zeros(numel(SensorsList),numel(SensorsList));
    SmoothnessGeodesic=zeros(numel(SensorsList),numel(SensorsList));
    SmoothnessAD=zeros(numel(SensorsList),numel(SensorsList));
    SmoothnessNCCA=zeros(numel(SensorsList),numel(SensorsList));
    SmoothnessKCCA=zeros(numel(SensorsList),numel(SensorsList));
    
    for SensorInd1=1:numel(SensorsList)
        K1=KernelsList1{SensorInd1};
        Datapoints1=DatapointsList1{SensorInd1};
        for SensorInd2=(1+SensorInd1):numel(SensorsList)
            K2=KernelsList2{SensorInd2};
            Datapoints2=DatapointsList2{SensorInd2};
            
            Dim=GetEffectiveDim(K1,K2,1e-3);
            [SNRVec,tVecCMR,t] = GetCMR(K1,K2,Dim,CMRParams);
            
            K_lt=(1-t)*K1+t*K2;
            K_gt=FixedGeodes( K1,K2,t,Dim );K_gt=real(K_gt);
            K_ad=0.5*( GetCS(K1)* GetCS(K2)'+ GetCS(K2)* GetCS(K1)');
%             [mU_ncca, mV_ncca] = NCCA(Datapoints1', Datapoints2', Lmax);
%             [mU_ncca, mV_ncca,Kx, Ky] = Kernel-NCCA(K1,K2,Datapoints1', Datapoints2', Lmax);
            [mU_ncca, mV_ncca,Kx, Ky] = NCCA_KernelsAsInputs(K1,K2, Lmax);
            if SmoothnessParams.UseGS
                [~,~,mU_kcca,mV_kcca] = KCCA(K1,K2,Lmax);
            else
                [mU_kcca,mV_kcca,~,~] = KCCA(K1,K2,Lmax);
            end
            
            s_l= GetSVecOfKernel(cX,K_lt,Lmax);
            s_g= GetSVecOfKernel(cX,K_gt,Lmax);
            s_ad= GetSVecOfKernel(cX,K_ad,Lmax);
            s_ncca=GetSVecsOfTwoEmbeddings(cX,mU_ncca,mV_ncca,Lmax);
            s_kcca=GetSVecsOfTwoEmbeddings(cX,mU_kcca,mV_kcca,Lmax);
            
            SmoothnessLinear(SensorInd1,SensorInd2)=s_l(SmoothnessParams.l);
            SmoothnessGeodesic(SensorInd1,SensorInd2)=s_g(SmoothnessParams.l);
            SmoothnessAD(SensorInd1,SensorInd2)=s_ad(SmoothnessParams.l);
            SmoothnessNCCA(SensorInd1,SensorInd2)=s_ncca(SmoothnessParams.l);
            SmoothnessKCCA(SensorInd1,SensorInd2)=s_kcca(SmoothnessParams.l);
        end
    end
    SmoothnessLinear=SmoothnessLinear+SmoothnessLinear';
    SmoothnessGeodesic=SmoothnessGeodesic+SmoothnessGeodesic';
    SmoothnessAD=SmoothnessAD+SmoothnessAD';
    SmoothnessNCCA=SmoothnessNCCA+SmoothnessNCCA';
    SmoothnessKCCA=SmoothnessKCCA+SmoothnessKCCA';
    
    diagonal=find(eye(size(SmoothnessLinear,1)));
    non_diagonal=find(1-(eye(size(SmoothnessLinear,1))));
    
    mean_geodesic=mean(SmoothnessGeodesic(non_diagonal));
    std_geodesic=std(SmoothnessGeodesic(non_diagonal));
    mean_linear=mean(SmoothnessLinear(non_diagonal));
    std_linear=std(SmoothnessLinear(non_diagonal));
    mean_ad=mean(SmoothnessAD(non_diagonal));
    std_ad=std(SmoothnessAD(non_diagonal));
    mean_ncca=mean(SmoothnessNCCA(non_diagonal));
    std_ncca=std(SmoothnessNCCA(non_diagonal));
    mean_kcca=mean(SmoothnessKCCA(non_diagonal));
    std_kcca=std(SmoothnessKCCA(non_diagonal));
    
    ScoresPerNoise=[ScoresPerNoise;[mean_geodesic,std_geodesic,mean_linear,std_linear,mean_ad,std_ad,mean_ncca,std_ncca,mean_kcca,std_kcca]];
end
fprintf([reverseStr]);

%% Print and save results
Methods=['Geodesic','Linear','AD','NCCA','KCCA'];
Noises=['Harmonic','Sawtooth'];
MetricsPerNoiseWithStd=[];
for row_ind=1:size(ScoresPerNoise,1)
    curr_row=[];
    for column_ind=1:2:size(ScoresPerNoise,2)
        curr_row=[curr_row,{sprintf('%.2f(%.2f)',ScoresPerNoise(row_ind,column_ind),ScoresPerNoise(row_ind,column_ind+1))}];
    end
    MetricsPerNoiseWithStd=[MetricsPerNoiseWithStd;curr_row];
end
test = array2table(MetricsPerNoiseWithStd,'VariableNames',{'Geodesic','Linear','AD','NCCA','KCCA'},'RowName',{'Harmonic','Sawtooth'});
disp(test)
writetable(test,fullfile(TablelsOutputFolder,"CM_Tabular_Results.txt"))


