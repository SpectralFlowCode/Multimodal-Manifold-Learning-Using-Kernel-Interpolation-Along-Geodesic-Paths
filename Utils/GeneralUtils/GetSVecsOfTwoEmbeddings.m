function [SVec] = GetSVecsOfTwoEmbeddings(cX,V,U,L)
        Lm=min(L,size(V,2));
        SVec1 = GetSVecsOfEmbeddings(cX,V,Lm);
        SVec2 = GetSVecsOfEmbeddings(cX,U,Lm);
        SVec=0.5*(SVec1+SVec2);
        if Lm<L
            SVec=[SVec,zeros(1,L-Lm)];
        end
end