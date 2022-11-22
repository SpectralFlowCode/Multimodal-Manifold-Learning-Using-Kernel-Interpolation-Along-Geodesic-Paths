function [ColumnStochasticK,K] = GetKernelFromData(Data,PreProcess,NormFac)
switch PreProcess
    case 'median'
        MedianBaseline=median(Data,1);Data=bsxfun(@minus,Data,MedianBaseline);
    case 'zscore'
        Data=zscore(Data);
    otherwise
end
tmpD=pdist2(Data',Data');
tmpA=exp(-tmpD/median(tmpD(:))/NormFac);
[ColumnStochasticK,K] =SingleIterationNorm(tmpA,1);
end

