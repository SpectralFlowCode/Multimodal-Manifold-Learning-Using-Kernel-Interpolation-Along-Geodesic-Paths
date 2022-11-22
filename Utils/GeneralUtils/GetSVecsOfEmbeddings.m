function [SVec] = GetSVecsOfEmbeddings(cX,V,L)
    V=V*diag(1./vecnorm(V, 2, 1));
    cXnorm=cX/norm(cX);
    L=min(L,size(V,2));
    SVec=cumsum((cXnorm*V(:,1:L)).^2);
end


