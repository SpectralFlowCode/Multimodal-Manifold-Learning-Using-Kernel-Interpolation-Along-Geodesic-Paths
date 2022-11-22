function [V,EigenVals] = GetSortedEVs(K,NumberOfEigenVals)
[V,D]=eigs(K,NumberOfEigenVals);EigenVals=diag(D);EigenVals=real(EigenVals);
[c,i]=sort(EigenVals,'descend');
EigenVals=EigenVals(i)';
V=real(V(:,i));

end

