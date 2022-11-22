disp(sprintf('\n--------- Generating figures for %s --------- ',FigPreamble))

%% Compute kernels
%projection using PCA (similarly as in the experiment in NCCA)
[~,projX1,~] = pca(X1');
[~,projX2,~] = pca(X2');

%Kernel contruction using nearest neighbors (similarly as in  the experiment in NCCA)
nn_kernel_params={};
nn_kernel_params.hx=2;
nn_kernel_params.nnx=20;
%
tmpA=GetKernelFromData_NearestNeighbors(projX1(:,1:100), 10,nn_kernel_params);
[K1,~] =SingleIterationNorm(tmpA,1);
%
tmpA=GetKernelFromData_NearestNeighbors(projX2(:,1:100), 10,nn_kernel_params);
[K2,~] =SingleIterationNorm(tmpA,1);

%% EVFD
Dim=GetEffectiveDim(K1,K2,EvfdParams.TolFac);
GeodesicKernelEigenValuesMat=[];%GeodesicKernelEigenValuesCorrMat=[];
GeodesicKernelEigenValuesCorrMat=[];
NumberOfEigenVals=10;

reverseStr = '';
tVec=linspace(0,1,EvfdParams.ntVec);
disp('Calculating the eigenvalues flow diagram');
for tind=1:length(tVec)
    percentDone = 100 * tind / length(tVec);
    msg = sprintf('   percentage done: %3.1f', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    t=tVec(tind);
    K_gt=FixedGeodes( K1,K2,t,Dim );K_gt=real(K_gt);
    %K_gt=(1-t)*K1+t*K2;
  
    [GeodesicV,EigenVals]=GetSortedEVs(K_gt,NumberOfEigenVals);
    GeodesicKernelEigenValuesMat=[GeodesicKernelEigenValuesMat;EigenVals];
end
fprintf([reverseStr]);
[GeodesicEigenValuesTupple,GeodesicCorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,GeodesicKernelEigenValuesMat,GeodesicKernelEigenValuesMat.^0,tVec);

%% Show EVFD
figure();
scatter((GeodesicEigenValuesTupple(:,1)),GeodesicEigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
sgtitle(sprintf('Eigenvalues flow diagram (EVFD)'),'FontSize', 40);
set(gcf,'Position',[489,180, 984,1110])
SaveFig(gcf,OutputFolder,FigPreamble+"_EVFD",1);

%% Show insets
Dim=GetEffectiveDim(K1,K2,EvfdParams.TolFac);
figure();

[V, ~, ~,~, ~,~] =DiffusionMapsFromKer( K1 , 1 );
figure()
scatter(V(:,2),V(:,3),50,n1,'filled')
axis off;axis tight;
title(sprintf('$v_{0}^{2}$ vs $v_{0}^{3}$'));hold on;
SaveFig(gcf,OutputFolder,FigPreamble+"_Inset_1",1);

[V, ~, ~,~, ~,~] =DiffusionMapsFromKer( K2 , 1 );
figure()
scatter(V(:,2),V(:,3),50,n2,'filled')
axis off;axis tight;
title(sprintf('$v_{1}^{2}$ vs $v_{1}^{3}$'));hold on;
SaveFig(gcf,OutputFolder,FigPreamble+"_Inset_0",1);


[SNRVec,tVecCMR,t] = GetCMR(K1,K2,Dim,CMRParams);
K=MnistFixedGeodes( K1,K2,t,Dim );K=real(K);
[V, ~, ~,singvals, ~,~] =DiffusionMapsFromKer( K , 1 );
figure()
scatter(V(:,2),V(:,3),50,labels,'filled')
axis off;axis tight;
title(sprintf('$v_{t^*}^{2}$ vs $v_{t^*}^{3}$'));hold on;
SaveFig(gcf,OutputFolder,FigPreamble+"_Inset_tstar",1);


