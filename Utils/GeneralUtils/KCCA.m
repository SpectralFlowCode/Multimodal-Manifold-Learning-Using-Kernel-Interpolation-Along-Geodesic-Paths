function [VV1,VV2,OVV1,OVV2] = KCCA(K1,K2,dim)
% VV1,VV2 - left and right kcca projections
% OVV1,OVV2 - left and right kcca projections after applying Gram-Schmidt

L   = size(K1, 1);
Z   = zeros(L);
kI  = 1e-3 * eye(L);
mG1 = [Z,K1*K2;
       K2*K1, Z];
mG1 = (mG1 + mG1') / 2;
mG2=[K1+kI,Z;
     Z, K2+kI];
[VV, LL]   = eigs(mG1, mG2, dim, 'largestreal');
[~, vIdx]  = sort(diag(LL), 'descend');
VV         = VV(:,vIdx);

VV  = real(VV);
VV1 = VV(1:L,:);
VV2 = VV(L+1:2*L,:);

% OVV1 = MyGramSchmidt(VV1);
% OVV2 = MyGramSchmidt(VV2);
OVV1 = GramSchmidt(VV1);
OVV2 = GramSchmidt(VV2);
%   figure(); imagesc(VV1'*VV1)
%   figure(); imagesc(OVV1'*OVV1)
end

