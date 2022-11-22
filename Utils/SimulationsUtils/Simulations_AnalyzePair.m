disp(sprintf('\n--------- Generating figures for %s --------- ',FigPreamble))

%% Compute kernels
D1 = pdist2(S1,S1,'euclidean');          
D2 = pdist2(S2,S2,'euclidean');           
eps1=median(D1(:))*KernelsParams.NormFac;
eps2=median(D2(:))*KernelsParams.NormFac;

switch KernelsParams.UseMutualScaleFac
    case 'Min'
        eps1=min([eps1,eps2]);eps2=eps1;
    case 'Mean'
        eps1=(eps1+eps2)/2;eps2=eps1;
end

A1=exp(-D1.^2/(eps1^2));
[ColumnStochasticK1,K1] =SingleIterationNorm(A1,1);
A2=exp(-D2.^2/(eps2^2));
[ColumnStochasticK2,K2] =SingleIterationNorm(A2,1);


%% Compute EVFDs 
Dim=GetEffectiveDim(K1,K2,2*TolFac);
LinearKernelEigenValuesMat=[];
GeodesicKernelEigenValuesMat=[];
reverseStr = '';
tVec=linspace(0,1,ntVec);
disp('Calculating the eigenvalues flow diagram');
for tind=1:length(tVec)
    percentDone = 100 * tind / length(tVec);
    msg = sprintf('   percentage done: %3.1f', percentDone);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    t=tVec(tind);
    K_gt=FixedGeodes( K1,K2,t,Dim );K_gt=real(K_gt);
%     K_gt=(1-t)*K1+t*K2;
    K_lt=(1-t)*K1+t*K2;
    
    [LinearV,EigenVals]=GetSortedEVs(K_lt,NumberOfEigenVals);
    LinearKernelEigenValuesMat=[LinearKernelEigenValuesMat;EigenVals];

    [GeodesicV,EigenVals]=GetSortedEVs(K_gt,NumberOfEigenVals);
    GeodesicKernelEigenValuesMat=[GeodesicKernelEigenValuesMat;EigenVals];
    
end
fprintf([reverseStr]);


%% Compare EVFDs (linear path vs geodesic path)
[LinearEigenValuesTupple,LinearCorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,LinearKernelEigenValuesMat,LinearKernelEigenValuesMat*0,tVec);
[GeodesicEigenValuesTupple,GeodesicCorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,GeodesicKernelEigenValuesMat,GeodesicKernelEigenValuesMat*0,tVec);

figure();
subplot(1,2,1);
scatter((GeodesicEigenValuesTupple(:,1)),GeodesicEigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
title('Geodesic($\gamma(t)$)','FontSize', 40);caxis([0 1]);colormap jet;
subplot(1,2,2);
scatter((LinearEigenValuesTupple(:,1)),LinearEigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
title('Linear($\mathbf{L}(t))$','FontSize', 40);caxis([0 1]);colormap jet;
sgtitle(sprintf('Eigenvalues flow diagram (EVFD)'),'FontSize', 40);
SaveFig(gcf,OutputFolder,FigPreamble+"_EVFDsComparison_SideBySide");

%% Mark analytical eigenvalues
CommonMus=0;
for i=0:20
        if strcmp(DataParams.Layout(1:6),'2DFlat')
            CommonMus=[CommonMus,(i*pi/2*DataParams.ScaleX)^2];
        end
        if strcmp(DataParams.Layout,'Torus_Polodial')
            CommonMus=[CommonMus,(floor(i/2)/DataParams.r1)^2];
%           CommonMus=[CommonMus,(i/DataParams.r1)^2];
        end
        if strcmp(DataParams.Layout,'Torus_Torodial')
            CommonMus=[CommonMus,(i/DataParams.R1)^2];
        end
end
CommonEigenvalsDiscrete1=exp(-CommonMus*(eps1^2)/4);


CommonMus=0;
for i=0:20
        if strcmp(DataParams.Layout(1:6),'2DFlat')
            CommonMus=[CommonMus,(i*pi/2*DataParams.ScaleX)^2];
        end
        if strcmp(DataParams.Layout,'Torus_Polodial')
           CommonMus=[CommonMus,(floor(i/2)/DataParams.r2)^2];
%             CommonMus=[CommonMus,(i/DataParams.r2)^2];
        end
        if strcmp(DataParams.Layout,'Torus_Torodial')
            CommonMus=[CommonMus,(i/DataParams.R2)^2];
        end
end
CommonEigenvalsDiscrete2=exp(-CommonMus*(eps2^2)/4);

%overlay on the same figure
figure();
scatter((GeodesicEigenValuesTupple(:,1)),GeodesicEigenValuesTupple(:,2),50,(0*GeodesicCorrEigenValuesRawStack),'b','filled','o','DisplayName','Eigenvalues flow diagram');
colormap jet;caxis([0 1]);%colorbar;
FontSize=35;xlFontSize=40;
set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', xlFontSize);
xlabel('$\log(\mu_t^i)$');ylabel('$t$')
hold on; scatter(log(CommonEigenvalsDiscrete1),CommonEigenvalsDiscrete1*0,500,'r','s','DisplayName',sprintf('Anal. eigenvalues ($\\mathcal{M}_x$)'));
hold on; scatter(log(CommonEigenvalsDiscrete2),ones(size(CommonEigenvalsDiscrete2)),500,'r','s','HandleVisibility','off');
for ind=1:20
    if ind==1
        h=line([log(CommonEigenvalsDiscrete1(ind)),log(CommonEigenvalsDiscrete2(ind))],[0,1],...
            'LineStyle','--','Color','r','Linewidth',5,'DisplayName',sprintf('Anal. trajectories ($\\mathcal{M}_x$)'));
    else
        h=line([log(CommonEigenvalsDiscrete1(ind)),log(CommonEigenvalsDiscrete2(ind))],[0,1],...
            'LineStyle','--','Color','r','Linewidth',5,'HandleVisibility','off');
    end
    h.Color(4) = 0.5;
end
xlim([min(GeodesicEigenValuesTupple(:,1)),max(GeodesicEigenValuesTupple(:,1))]);
hLegend = findobj(gcf, 'Type', 'Legend');
lgd=legend('Location','eastoutside');
lgd.FontSize=15;
SaveFig(gcf,OutputFolder,FigPreamble+"_EVFD_WithAnalyticEigenvalues");


%% Show diffusion patterns
% Compute kernels for visualizations
eps1=median(D1(:))*KernelsParams.NormFacVisualization;
eps2=median(D2(:))*KernelsParams.NormFacVisualization;
A1=exp(-D1.^2/(eps1^2));
[ColumnStochasticK1,K1] =SingleIterationNorm(A1,1);
A2=exp(-D2.^2/(eps2^2));
[ColumnStochasticK2,K2] =SingleIterationNorm(A2,1);


if strcmp(DataParams.Layout(1:6),'2DFlat')
    maxL=1.2*max([DataParams.ScaleX,DataParams.ScaleY,DataParams.ScaleZ]);
    PointSize=10;

else
    max1=1*(DataParams.r1+DataParams.R1);
    max2=1*(DataParams.r2+DataParams.R2);
    PointSize=100;

end

NumberOfSteps=1;
v0=[1,zeros(1,N-1)];
tVecToShow=[0,0.2,0.5,1];
Dim=GetEffectiveDim(K1,K2,TolFac);
rot=[  -68.9000   67.2980];

figure();
ind=1;
for tToShow=tVecToShow
    K=FixedGeodes( K1,K2,tToShow,Dim);K=real(K);%figure(); imagesc(K(I,I));
    %K=(1-tToShow)*K1+tToShow*K2;
    vt=real(K^NumberOfSteps*v0');   
    if strcmp(DataParams.Layout(1:6),'2DFlat')
        subplot(4,2,ind);
        scatter(S1(:,1),S1(:,2),PointSize,vt,'filled');axis(maxL*[-1,1,-1,1]);
        axis off;     
        subplot(4,2,ind+1);
        scatter(S2(:,1),S2(:,2),PointSize,vt,'filled');axis(maxL*[-1,1,-1,1]);
        axis off;   
    else
        subplot(4,2,ind);
        scatter3(S1(:,1),S1(:,2),S1(:,3),PointSize,vt,'filled');axis(max1*[-1,1,-1,1]);
        axis off;     axis tight;    view(rot);
        subplot(4,2,ind+1);
        scatter3(S2(:,1),S2(:,2),S2(:,3),PointSize,vt,'filled');axis(max2*[-1,1,-1,1]);
        axis off;      axis tight;    view(rot);
    end
    ind=ind+2;
end
set(gcf,'Position',[343,198,717,966]);
SaveFig(gcf,OutputFolder,FigPreamble+"_DiffusionPatterns",-1);



