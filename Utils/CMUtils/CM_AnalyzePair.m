disp(sprintf('\n--------- Generating figures for %s --------- ',FigPreamble))

%% Load samples
Target=ConditionNames{TargetInd};
%filter outliers
if TargetInd==1
    %when predicting the cooler, in order to preserve the diversity of the
    %labels we preserve the beginnig of the samples, and filter outliers
    Data=Samples('EPS1');Data=Data(:,RelevantInds);
    MedianBaseline=median(Data,1);Data=bsxfun(@minus,Data,MedianBaseline);
    tmpD=pdist2(Data',Data');
    e=abs(tmpD(1,:)-median(tmpD(1,:)));
    medstd=median(e);
    %figure(); plot(e);hold on; line([1,size(Data,2)],3*[medstd,medstd]);
    outliers=find(e>3*medstd);
    IndsToProcess=setdiff(RelevantInds,RelevantInds(outliers));
else
    %remove outliers from the beginning of the sample
    IndsToProcess=RelevantInds(460:end);
end
GT=Conditions(Target);GT=GT(IndsToProcess);

%% Load kernels
[csK1,K1] = GetKernelFromNoisySensor(Samples,SensorsToLoad{SensorInd1},PreProcess,IndsToProcess,NormFac1,noise1,NoiseParams1,false);
[csK2,K2] = GetKernelFromNoisySensor(Samples,SensorsToLoad{SensorInd2},PreProcess,IndsToProcess,NormFac2,noise2,NoiseParams2,false);

Sensor1=SensorsToLoad{SensorInd1};
Sensor2=SensorsToLoad{SensorInd2};

%% Compute EVFD + smoothness evaluation metrics
Lmax=100;
cX=GT-mean(GT);
DX= pdist2(GT',GT','euclidean');

Dim=GetEffectiveDim(K1,K2,th);
LinearKernelEigenValuesMat=[];%LinearKernelEigenValuesCorrMat=[];
GeodesicKernelEigenValuesMat=[];%GeodesicKernelEigenValuesCorrMat=[];
LinearKernelEigenValuesCorrMat=[];
GeodesicKernelEigenValuesCorrMat=[];
reverseStr = '';
tVec=linspace(0,1,ntVec);
disp('Calculating the eigenvalues flow diagram');
for tind=1:length(tVec)
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

%% Compute smoothness evaluation metrics for AD
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
[SNRVec,tVecCMR,tstar] = GetCMR(K1,K2,Dim,figCMRParams);
inset1=2;inset2=3;inset3=4;
if strcmp(FigPreamble,'Figure2')
    for t=[0,tstar,1]
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
        title(sprintf('Top2 dominant comoponents\n at $t=%s$ ',t));hold on;
        SaveFig(gcf,OutputFolder,FigPreamble+"_Inset_"+t(1));
    end
else
    t_insets=[0,tstar,1];
    CameraPositions=[ [4.8869,0.4469,7.6869];[-0.4488,-5.7611,2.3454];[-0.4488,-5.7611,2.3454]];
    for t_ind=1:3
        t=t_insets(t_ind);
        K=FixedGeodes( K2,K1,t,Dim );K=real(K);
        [V, ~, ~,singvals, ~,~] =DiffusionMapsFromKer( K , 1 );
        figure();
        scatter3(V(:,inset1),V(:,inset2),V(:,inset3),50,GT,'filled')
        axis off;axis tight;
        if not(t==1 | t==0)
            t='t^*';
        else
            t=num2str(t);
        end
        title(sprintf('Top3 dominant comoponents\n at $t=%s$ ',t));hold on;
        set(gca, 'CameraPosition',CameraPositions(t_ind,:));
        SaveFig(gcf,OutputFolder,FigPreamble+"_Inset_"+t(1));
    end
end
    