function [ColumnStochasticK,K] =SingleIterationNorm(W,alpha)

A=diag(1./(sum(W).^alpha))*W*diag(1./(sum(W).^alpha));
K=diag((1./sum(A)).^(1/2))*A*diag((1./sum(A)).^(1/2));
ColumnStochasticK=diag((1./sum(A)).^(-1/2))*K*diag((1./sum(A)).^(1/2));


end


