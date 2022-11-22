% clc;
% clear all;
% close all;
addpath(genpath('Utils'));

%% Load data and parameters
ENose_GetParameters;
ENose_GetData;

%% Evaluate all pairs
disp(sprintf('Generating E-Nose objective results'))
l=SmoothnessParams.l;
Lmax=SmoothnessParams.l;
ScoresPerGas=[];
reverseStr = '';
for TargetInd=1:6
    msg = sprintf('Evaluating %s (%g/6)', GasNames{TargetInd},TargetInd);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    SensorsList=DatasetsParams.SensorsListPerGas{TargetInd};
    Target=GasNames{TargetInd};
    GT=GasGT(Target)';GT=GT(RelevantInds);
    cX=GT-mean(GT);
    Samples=SamplesCell{TargetInd};
    
    KernelsList={};DatapointsList={};
    for SensorInd=SensorsList
        SensorToUse=SensorsToLoad{SensorInd};
        [ColumnStochasticK,K,Data] = GetKernelFromSensor(Samples,SensorToUse,SmoothnessParams.PreProcess,RelevantInds,SmoothnessParams.NormFac);
        KernelsList=[KernelsList,K];
        DatapointsList=[DatapointsList,Data];
    end
    %compute smoothmess for each pair
    SVecsMatLinear=zeros(numel(SensorsList),numel(SensorsList),Lmax);
    SVecsMatGeodesic=zeros(numel(SensorsList),numel(SensorsList),Lmax);
    SVecsMatAD=zeros(numel(SensorsList),numel(SensorsList),Lmax);
    SVecsMatNCCA=zeros(numel(SensorsList),numel(SensorsList),Lmax);
    SVecsMatKCCA=zeros(numel(SensorsList),numel(SensorsList),Lmax);
    
    for SensorInd1=1:numel(SensorsList)
        K1=KernelsList{SensorInd1};
        Datapoints1=DatapointsList{SensorInd1};
        
        for SensorInd2=1:numel(SensorsList)
            K2=KernelsList{SensorInd2};
            Datapoints2=DatapointsList{SensorInd2};
            Dim=GetEffectiveDim(K1,K2,1e-3);
            [SNRVec,tVecCMR,t] = GetCMR(K1,K2,Dim,CMRParams);
            
            K_lt=(1-t)*K1+t*K2;
            K_gt=FixedGeodes( K1,K2,t,Dim );K_gt=real(K_gt);
            K_ad=0.5*( GetCS(K1)* GetCS(K2)'+ GetCS(K2)* GetCS(K1)');
            [mU_ncca, mV_ncca,Kx, Ky] = NCCA_KernelsAsInputs(K1,K2, Lmax);

            if SmoothnessParams.UseGS
                [~,~,mU_kcca,mV_kcca] = KCCA(K1,K2,Lmax);
            else
                [mU_kcca,mV_kcca,~,~] = KCCA(K1,K2,Lmax);
            end
            
            SVecsMatLinear(SensorInd1,SensorInd2,:)= GetSVecOfKernel(cX,K_lt,Lmax);
            SVecsMatGeodesic(SensorInd1,SensorInd2,:)= GetSVecOfKernel(cX,K_gt,Lmax);
            SVecsMatAD(SensorInd1,SensorInd2,:)=GetSVecOfKernel(cX,K_ad,Lmax);
            SVecsMatNCCA(SensorInd1,SensorInd2,:)=GetSVecsOfTwoEmbeddings(cX,mU_ncca,mV_ncca,Lmax);
            SVecsMatKCCA(SensorInd1,SensorInd2,:)=GetSVecsOfTwoEmbeddings(cX,mU_kcca,mV_kcca,Lmax);     
        end
    end
    
    [tmp_l,tmp_g,tmp_ad,tmp_ncca,tmp_kcca]=...
        deal(SVecsMatLinear(:,:,l),SVecsMatAD(:,:,l),SVecsMatGeodesic(:,:,l),SVecsMatNCCA(:,:,l),SVecsMatKCCA(:,:,l));
    diagonal=find(eye(size(tmp_g,1)));
    non_diagonal=find(1-(eye(size(tmp_g,1))));
    %
    mean_geodesic=mean(tmp_g(non_diagonal));
    mean_linear=mean(tmp_l(non_diagonal));
    mean_ad=mean(tmp_ad(non_diagonal));
    mean_ncca=mean(tmp_ncca(non_diagonal));
    mean_kcca=mean(tmp_kcca(non_diagonal));
    %
    std_geodesic=std(tmp_g(non_diagonal));
    std_linear=std(tmp_l(non_diagonal));
    std_ad=std(tmp_ad(non_diagonal));
    std_ncca=std(tmp_ncca(non_diagonal));
    std_kcca=std(tmp_kcca(non_diagonal));
    
    ScoresPerGas=[ScoresPerGas;[mean_geodesic,std_geodesic,mean_linear,std_linear,mean_ad,std_ad,mean_ncca,std_ncca,mean_kcca,std_kcca]];
end
fprintf([reverseStr]);

%% Print and save results
ReorderedMetricsPerGas=ScoresPerGas([3,4,5,1,2,6],:);
ScoresWithStd=[];
for row_ind=1:size(ScoresPerGas,1)
    curr_row=[];
    for column_ind=1:2:size(ReorderedMetricsPerGas,2)
        curr_row=[curr_row,{sprintf('%.2f(%.2f)',ReorderedMetricsPerGas(row_ind,column_ind),ReorderedMetricsPerGas(row_ind,column_ind+1))}];
    end
    ScoresWithStd=[ScoresWithStd;curr_row];
end
mean_row=[];
for column_ind=1:2:size(ScoresPerGas,2)
    mean_row=[mean_row,{sprintf('%.2f(%.2f)',mean(ScoresPerGas(:,column_ind)),mean(ScoresPerGas(:,1+column_ind)))}];
end
ScoresWithStd=[ScoresWithStd;mean_row];


test = array2table(ScoresWithStd,'VariableNames',{'Geodesic','Linear','AD','NCCA','KCCA'},'RowName',[GasNames([3,4,5,1,2,6]),'Mean']);
disp(test)
writetable(test,fullfile(TablelsOutputFolder,"ENose_Tabular_Results.txt"))




