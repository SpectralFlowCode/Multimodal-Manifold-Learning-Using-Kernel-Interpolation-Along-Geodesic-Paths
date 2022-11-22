disp(sprintf('\n--------- Generating figures for %s --------- ',FigPreamble))

%% Load samples
Target=GasNames{TargetInd};
GT=GasGT(Target)';GT=GT(RelevantInds);
Samples=SamplesCell{TargetInd};

[csK1,K1] = GetKernelFromSensor(Samples,SensorsToLoad{SensorInd1},PreProcess,RelevantInds,NormFac1);
[csK2,K2] = GetKernelFromSensor(Samples,SensorsToLoad{SensorInd2},PreProcess,RelevantInds,NormFac2);
Sensor1=SensorsToLoad{SensorInd1};
Sensor2=SensorsToLoad{SensorInd2};

%% EVFD + smoothness evaluation metrics (as function of t)
Lmax=100;
cX=GT-mean(GT);
DX= pdist2(GT',GT','euclidean');

Dim=GetEffectiveDim(K1,K2,1e-3);


LinearKernelEigenValuesMat=[];%LinearKernelEigenValuesCorrMat=[];
GeodesicKernelEigenValuesMat=[];%GeodesicKernelEigenValuesCorrMat=[];
LinearKernelEigenValuesCorrMat=[];
GeodesicKernelEigenValuesCorrMat=[];

reverseStr = '';
tVec=linspace(0,1,ntVec);
disp('Calculating the eigenvalues flow diagram');
for tind=1:length(tVec)
    %     tind
    percentDone = 100 * tind / length(tVec);
    msg = sprintf('   percentage done: %3.1f', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    t=tVec(tind);
    K_gt=FixedGeodes( K1,K2,t,Dim );K_gt=real(K_gt);
    K_lt=(1-t)*K1+t*K2;
    
    [LinearV,EigenVals]=GetSortedEVs(K_lt,NumberOfEigenVals);
    LinearKernelEigenValuesMat=[LinearKernelEigenValuesMat;EigenVals];
    LinearKernelEigenValuesCorrMat=[LinearKernelEigenValuesCorrMat;abs(corr(GT',LinearV(:,1:NumberOfEigenVals)))];
    
    [GeodesicV,EigenVals]=GetSortedEVs(K_gt,NumberOfEigenVals);
    GeodesicKernelEigenValuesMat=[GeodesicKernelEigenValuesMat;EigenVals];
    GeodesicKernelEigenValuesCorrMat=[GeodesicKernelEigenValuesCorrMat;abs(corr(GT',GeodesicV(:,1:NumberOfEigenVals)))];
end
fprintf([reverseStr]);

%%  Compute smoothness evaluation metrics for AD
K_ad=0.5*( GetCS(K1)* GetCS(K2)'+ GetCS(K2)* GetCS(K1)'); %BackwardForward implementation

%% Show EVFDs
[LinearEigenValuesTupple,LinearCorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,LinearKernelEigenValuesMat,LinearKernelEigenValuesCorrMat,tVec);
[GeodesicEigenValuesTupple,GeodesicCorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,GeodesicKernelEigenValuesMat,GeodesicKernelEigenValuesCorrMat,tVec);


figure();
subplot(1,2,1);
scatter((GeodesicEigenValuesTupple(:,1)),GeodesicEigenValuesTupple(:,2),50,GeodesicCorrEigenValuesRawStack,'filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
title('Geodesic($\gamma(t)$)','FontSize', 40);caxis([0 1]);colormap jet;
subplot(1,2,2);
scatter((LinearEigenValuesTupple(:,1)),LinearEigenValuesTupple(:,2),50,LinearCorrEigenValuesRawStack,'filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
title('Linear($\mathbf{L}(t))$','FontSize', 40);caxis([0 1]);colormap jet;
sgtitle(sprintf('Evfds connecting sensor %s and sensor %s\n colored by %s',Sensor1,Sensor2,Target),'FontSize', 40);
SaveFig(gcf,OutputFolder,FigPreamble+"_EVFD");


%% Show insets
for curr_inset=common_insets
    curr_inset=curr_inset{:};
    [t,inset1,inset2]=deal(curr_inset(1),curr_inset(2),curr_inset(3));
    K=FixedGeodes( K1,K2,t,Dim );K=real(K);
    [V, ~, ~,singvals, ~,~] =DiffusionMapsFromKer( K , 1 );
    figure();
    scatter(V(:,inset1),V(:,inset2),50,GT,'filled')
    axis off;axis tight;
    if not(t==1 | t==0)
        t='t^*';
    else
        t=num2str(t);
    end
    title(sprintf('Common components: \n $v_{%s}^{%g}$ vs $v_{%s}^{%g}$',t,inset1,t,inset2));hold on;
    SaveFig(gcf,OutputFolder,FigPreamble+"_CommonInset_"+t(1));
end

for curr_inset=noise_insets
    curr_inset=curr_inset{:};
    [t,inset1,inset2]=deal(curr_inset(1),curr_inset(2),curr_inset(3));
    K=FixedGeodes( K1,K2,t,Dim );K=real(K);
    [V, ~, ~,singvals, ~,~] =DiffusionMapsFromKer( K , 1 );
    figure();
    scatter(V(:,inset1),V(:,inset2),50,GT,'filled')
    axis off;axis tight;
    title(sprintf('Measurement-specific components: \n $v_{%g}^{%g}$ vs $v_{%g}^{%g}$',t,inset1,t,inset2));hold on;
        SaveFig(gcf,OutputFolder,FigPreamble+"_NoiseInset_"+num2str(t));
end


