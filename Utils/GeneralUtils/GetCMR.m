function [SNRVec,tVec,tstar] = GetCMR(K1,K2,Dim,Params)
v2struct(Params)
tVec=linspace(0,1,ntVec);

%--- calculate using multiplication along the path

% GlobalProdK=eye(size(K1,1));
% for t=linspace(0.3,0.7,Np)
%     if Fast
%         K=(1-t)*K1+t*K2;
%     else
%         K=FixedGeodes( K1,K2,t,Dim );K=real(K);
%     end
%     GlobalProdK=GlobalProdK*K;
% end
% [VG,EigenValsG]=GetSortedEVs(GlobalProdK,n);

%--- calculate using average smoothness maximization along the path
%(according to: Dietrich, Felix, et al. "Spectral discovery of jointly smooth features for multimodal data." 
% SIAM Journal on Mathematics of Data Science 4.1 (2022): 410-430.

VGs=[];
for t=acc_path
    if Fast
        K=(1-t)*K1+t*K2;
    else
        K=FixedGeodes( K1,K2,t,Dim );K=real(K);
    end
    [VG,~]=GetSortedEVs(K,n);
    VGs=[VGs,VG];
end
% [coeff,score,latent] = pca(VGs);
% VG=score(:,1:n);
[U,S,V]=svd(VGs*VGs');
VG=U(:,1:n);


SNRVec=[];
for CurrT=linspace(0,1,ntVec)
    %CurrT
    if Fast
        currK=(1-CurrT)*K1+CurrT*K2;
    else
        currK=FixedGeodes( K1,K2,CurrT,Dim );currK=real(currK);
    end
    [V,EigenVals]=GetSortedEVs(currK,Ne);
    
    SpreadingMat=abs(V'*VG);
    % figure(); imagesc(SpreadingMat(2:end,2:end))
    % figure(); myplot(SpreadingMat(:,2:end))
    
    
    ProbVec=abs(SpreadingMat(2:end,2:n));
    ProbVec=diag(1./sum(ProbVec,2))*ProbVec;
    % figure(); imagesc(ProbVec)
    EntVec=zeros(1,Ne-1);
    MaxMinVec=zeros(1,Ne-1);
    for ind=1:Ne-1
        tmp=nonzeros(ProbVec(ind,:));
        EntVec(ind)=-sum(tmp.*log(tmp));
        MaxMinVec(ind)=max(tmp)-min(tmp);
    end
    
    Scores=(log(n-1)-EntVec);
    [c,SortedInds]=sort(Scores,'descend');
    SignalInds=SortedInds(1);
    NoiseInds=setdiff(SortedInds,SignalInds);
    EntTopKCMR=sum(EigenVals(1+SignalInds))/sum(EigenVals(1+NoiseInds));
    
    NormalizedScores=Scores/sum(Scores);
    tmp_wsignal=sum(EigenVals(2:Ne).*NormalizedScores);
    tmp_wnoise=sum(EigenVals(2:Ne).*(1-NormalizedScores));
    EntWeightedCMR=tmp_wsignal/tmp_wnoise;
    
    tmp_wsignal=sum(EigenVals(2:Ne).*(log(n-1)-EntVec));
    tmp_wnoise=sum(EigenVals(2:Ne).*(EntVec));
    EntCMR=tmp_wsignal/tmp_wnoise;
    
    
    Scores=MaxMinVec;
    [c,SortedInds]=sort(Scores,'descend');
    SignalInds=SortedInds(1);
    NoiseInds=setdiff(SortedInds,SignalInds);
    nScores=Scores(1:end);
    nScores=nScores/sqrt(sum(nScores.^2));
    nE=EigenVals(2:(Ne));
    nE=nE/sqrt(sum(nE.^2));
    tmp_wsignal=sum(nE.*nScores');
    tmp_wnoise=sum(nE.*(1-nScores'));
    MaxMinCMR=tmp_wsignal/tmp_wnoise;
    
    Smoothness=trace(V'*VG*VG'*V)/trace(V'*V); 
    Scores=diag(V(:,2:end)'*VG*VG'*V(:,2:end))';
    NormalizedScores=Scores/sum(Scores);
    tmp_wsignal=sum(EigenVals(2:Ne).*NormalizedScores);
    tmp_wnoise=sum(EigenVals(2:Ne).*(1-NormalizedScores));
    WeightedSmoothness=tmp_wsignal/tmp_wnoise;

    SNRVec=[SNRVec,[EntTopKCMR;EntWeightedCMR;EntCMR;MaxMinCMR;Smoothness;WeightedSmoothness]];
end


switch criterion
    case 'ent_topk'
        criterion_ind=1;
    case 'ent_weighted'
        criterion_ind=2;
    case 'ent'
        criterion_ind=3;
    case 'max_min'
        criterion_ind=4;
    case 'smoothness'
        criterion_ind=5;
    case 'smoothness_weighted'
        criterion_ind=6;
end
[~,i]=max(SNRVec(criterion_ind,:));
tstar=tVec(i);
% figure(); myplot(zscore(SNRVec'));
% figure(); plot(SNRVec(5,:));


end

