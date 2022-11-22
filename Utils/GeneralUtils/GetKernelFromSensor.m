function [ColumnStochasticK,K,Data] = GetKernelFromSensor(Samples,SensorToUse,PreProcess,RelevantInds,NormFac)
Data=Samples(SensorToUse);Data=Data(:,RelevantInds);
switch PreProcess
    case 'none'
    case 'median'
        MedianBaseline=median(Data,1);
        Data=bsxfun(@minus,Data,MedianBaseline);
    case 'zscore'
        Data=zscore(Data')';
%         figure(); myplot(zscore(Data')); %if zscored
%         figure(); myplot(Data'); %if not
    otherwise
        print('Bug')
end
tmpD=pdist2(Data',Data');
tmpA=exp(-tmpD/median(tmpD(:))/NormFac);
[ColumnStochasticK,K] =SingleIterationNorm(tmpA,1);
end

