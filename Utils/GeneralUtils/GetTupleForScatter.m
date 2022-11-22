function [LinearEigenValuesTupple,LinearCorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,LinearKernelEigenValuesMat,LinearKernelEigenValuesCorrMat,tVec)
tmpLinear=log(LinearKernelEigenValuesMat(:,2:NumberOfEigenVals)');
tmpLinearCorr=(LinearKernelEigenValuesCorrMat(:,2:NumberOfEigenVals)');
LinearEigenValuesRawStack=tmpLinear(:);
LinearCorrEigenValuesRawStack=tmpLinearCorr(:);
tRawStack=kron(tVec(1:size(tmpLinear,2)),ones(1,size(tmpLinear,1)))';
LinearEigenValuesTupple=[LinearEigenValuesRawStack,tRawStack];

end

