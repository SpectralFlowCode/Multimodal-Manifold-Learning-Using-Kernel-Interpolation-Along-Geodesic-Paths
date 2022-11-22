function csK = GetCS(K)
%returns column stochastic matrix
D=sum(K,1);
csK=K*diag(D.^(-1));
end

