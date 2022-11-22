
function [SVec] = GetSVecOfKernel(cX,K,L)
        [V,EigenVals]=GetSortedEVs(K,L);
        cXnorm=cX/norm(cX);
        SVec=cumsum((cXnorm*V(:,1:L)).^2);   
end