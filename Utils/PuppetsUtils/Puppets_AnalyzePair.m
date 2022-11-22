disp(sprintf('\n--------- Generating figures for %s --------- ',FigPreamble))

%% Compute kernels
d1 = pdist2(s1,s1,'euclidean');             % The distance between samples, based only on snapshots taken by Camera 1.
d2 = pdist2(s2,s2,'euclidean');             % The distance between samples, based only on snapshots taken by Camera 2.
D1=sqrt(d1);
D2=sqrt(d2);

UseMutualEps='none';
eps1=median(D1(:))*NormFac;
eps2=median(D2(:))*NormFac;

switch UseMutualEps
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
GT=MapEmbdbulldog(:,2)';
cX=GT-mean(GT);

Lmax=100;
PolyNormFac=Inf;
Dim=GetEffectiveDim(K1,K2,TolFac);
tVec=linspace(0,1,ntVec);

LinearKernelEigenValuesMat=[];
GeodesicKernelEigenValuesMat=[];
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

%% AD
K_ad=0.5*( GetCS(K1)* GetCS(K2)+ GetCS(K2)* GetCS(K1));%K=K^2;

%% Show EVFDs side by side
[LinearEigenValuesTupple,LinearCorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,LinearKernelEigenValuesMat,LinearKernelEigenValuesCorrMat,tVec);
[GeodesicEigenValuesTupple,GeodesicCorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,GeodesicKernelEigenValuesMat,GeodesicKernelEigenValuesCorrMat,tVec);

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

%% Overlay EVFDs
figure();
scatter((GeodesicEigenValuesTupple(:,1)),GeodesicEigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
hold on;
scatter((LinearEigenValuesTupple(:,1)),LinearEigenValuesTupple(:,2),50,'r','filled','o','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
legend('$\gamma(t)$','$\mathbf{L}(t))$');
sgtitle(sprintf('A comparison between the EVFD obtained using a geodesic  \n interpolation ($\\gamma(t)$) and the EVFD obtained using linear interpolation ($\\mathbf{L}(t))$)'),'FontSize', 30);
SaveFig(gcf,OutputFolder,FigPreamble+"_EVFDsComparison_Overlay");



%% Show insets
FontSize=10;xlFontSize=60;ylFontSize=60;
PointSize=100;

for t0=[0,0.5,1]
    K=FixedGeodes( K1,K2,t0,Dim );K=real(K);
    [MapEmbd, ~, ~,singvals, ~,~] =DiffusionMapsFromKer( K , 1);
    
    figure();
    subplot(3,3,1);
    scatter(AngleYoda,MapEmbd(:,2),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{2}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{2}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none');
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,2);
    scatter(AngleBulldog,MapEmbd(:,2),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{2}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{2}_{t^{*}}$',t0));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,3);
    scatter(AngleBunny,MapEmbd(:,2),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{2}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{2}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    
    %
    subplot(3,3,4);
    scatter(AngleYoda,MapEmbd(:,3),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{3}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{2}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,5);
    scatter(AngleBulldog,MapEmbd(:,3),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{3}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{3}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,6);
    scatter(AngleBunny,MapEmbd(:,3),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{3}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{3}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    %
    subplot(3,3,7);
    scatter(AngleYoda,MapEmbd(:,4),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{4}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{4}_{t^{*}}$'));
    end
    %     if t0==0
    xlabel('$\theta_{\mathrm{yoda}}$');
    %     end
    ax.YLabel.Visible = 'on';
    ax.XLabel.Visible = 'on';
    
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    axis off;
    
    subplot(3,3,8);
    scatter(AngleBulldog,MapEmbd(:,4),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{4}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{4}_{t^{*}}$',t0));
    end
    xlabel('$\theta_{\mathrm{bulldog}}$');
    ax.YLabel.Visible = 'on';
    ax.XLabel.Visible = 'on';
    
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,9);
    scatter(AngleBunny,MapEmbd(:,4),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{4}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{4}_{t^{*}}$'));
    end
    xlabel('$\theta_{\mathrm{bunny}}$');
    ax.YLabel.Visible = 'on';
    ax.XLabel.Visible = 'on';
    
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    SaveFig(gcf,OutputFolder,[FigPreamble,sprintf('_Inset_%g',t0)]);
end