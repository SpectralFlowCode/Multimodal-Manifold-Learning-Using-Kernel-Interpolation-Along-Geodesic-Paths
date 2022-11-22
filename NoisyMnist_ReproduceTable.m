% clc;
% clear all;
% close all;
addpath(genpath('Utils'));

%% Load data and parameters
NoisyMnist_GetParameters;
disp(sprintf('Generating tables for %s',FigPreamble))


%% Initializations
nccas=[];
kccas=[];
ours=[];
linears=[];
ads=[];

%% Perform simulations
for ind=1:SimParams.nReps
    %% Load data and partition to train-test
    NoisyMnist_GetData;
    cv = cvpartition(N,'HoldOut',0.25);
    train_inds=find(~cv.test);
    test_inds=find(cv.test);
    
    %% Compute kernels
    [~,projX1,~] = pca(X1');
    [~,projX2,~] = pca(X2');
    
    tmpD=pdist2(projX1(:,1:100),projX1(:,1:100)).^2;
    tmpA=exp(-tmpD/median(tmpD(:))/NormFac);
    [ColumnStochasticK,K1] =SingleIterationNorm(tmpA,1);
    
    tmpD=pdist2(projX2(:,1:100),projX2(:,1:100)).^2;
    tmpA=exp(-tmpD/median(tmpD(:))/NormFac);
    [ColumnStochasticK,K2] =SingleIterationNorm(tmpA,1);
    
   
    Dim=GetEffectiveDim(K1,K2,EvfdParams.TolFac);     
    [SNRVec,tVecCMR,t] = GetCMR(K1,K2,Dim,CMRParams);
%     t=0.5;
    
    %% Using geodesic
    K=FixedGeodes( K1,K2,t,Dim );K=real(K);
    [V,EigenVals]=GetSortedEVs(K,EvfdParams.NumberOfEigenVals);
    Model_ours = fitcecoc(V(train_inds,:),labels(train_inds)','Learners',templateSVM('Standardize',true));
    [label,score] = predict(Model_ours,V(test_inds,:));
    acc_ours=sum((label-labels(test_inds)==0))/numel(test_inds);
    
    %% Using linear
    K=(1-t)*K1+t*K2;
    [V,EigenVals]=GetSortedEVs(K,EvfdParams.NumberOfEigenVals);
    Model_linear = fitcecoc(V(train_inds,:),labels(train_inds)','Learners',templateSVM('Standardize',true));
    [label,score] = predict(Model_linear,V(test_inds,:));
    acc_linear=sum((label-labels(test_inds)==0))/numel(test_inds);
    
    
    %% Using KCCA
    [mU_kcca,mV_kcca,~,~] = KCCA(K1,K2,EvfdParams.NumberOfEigenVals);
    Model_kcca = fitcecoc(mU_kcca(train_inds,:),labels(train_inds)','Learners',templateSVM('Standardize',true));
    [label,score] = predict(Model_kcca,mU_kcca(test_inds,:));
    acc_kcca=sum((label-labels(test_inds))==0)/numel(test_inds);

    %% Using NCCA
    [mU_ncca, mV_ncca,Kx, Ky] = NCCA_KernelsAsInputs(K1,K2,  EvfdParams.NumberOfEigenVals);
    Model_ncca = fitcecoc(mU_ncca(train_inds,:),labels(train_inds)','Learners',templateSVM('Standardize',true));
    [label,~] = predict(Model_ncca,mU_ncca(test_inds,:));
    acc_ncca=sum((label-labels(test_inds))==0)/numel(test_inds);
    
    %% Using AD
    K_ad=0.5*( GetCS(K1)* GetCS(K2)'+ GetCS(K2)* GetCS(K1)');
    [V,EigenVals]=GetSortedEVs(K_ad, EvfdParams.NumberOfEigenVals);
    Model_ad = fitcecoc(V(train_inds,:),labels(train_inds)','Learners',templateSVM('Standardize',true));
    [label,score] = predict(Model_ad,V(test_inds,:));
    acc_ad=sum((label-labels(test_inds)==0))/numel(test_inds);
    
    %% Concatenate results
    nccas=[nccas,acc_ncca];
    kccas=[kccas,acc_kcca];
    ours=[ours,acc_ours];
    linears=[linears,acc_linear];
    ads=[ads,acc_ad];
end

%% Aggregate results
mean_nccas=mean(nccas);
mean_kccas=mean(kccas);
mean_ours=mean(ours);
mean_linears=mean(linears);
mean_ads=mean(ads);
Table=[{ours,linears,ads,nccas,kccas}];

%% Compute STDs and print aggregaed results
AggTableText=cell(size(Table,1),size(Table,2));
AggTable=zeros(size(Table,1),2*size(Table,2));
for row_ind=1:size(Table,1)
    text=[''];
    for col_ind=1:size(Table,2)
        AggTable(row_ind,2*(col_ind-1)+1)=mean(Table{row_ind,col_ind});
        AggTable(row_ind,2*(col_ind))=std(Table{row_ind,col_ind});
        AggTableText{row_ind,col_ind}=sprintf('$%.2f (%.2f)$',mean(Table{row_ind,col_ind}),std(Table{row_ind,col_ind}));
    end
    text=[text,'\\\\'];
end
sympref('FloatingPointOutput',1);
latex_table = latex(sym(AggTable));
tmp=array2table(AggTableText(:,1:end),'VariableNames',{'Geodesic','Linear','AD','NCCA','KCCA'});
disp(tmp);
writetable(tmp,fullfile(TablelsOutputFolder,FigPreamble+"_Tabular_Results.txt"))

