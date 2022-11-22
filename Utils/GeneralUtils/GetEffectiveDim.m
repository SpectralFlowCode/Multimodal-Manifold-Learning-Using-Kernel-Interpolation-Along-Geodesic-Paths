function [Dim] = GetEffectiveDim(K1,K2,th)
    [~,D1]=eig(K1);[~,D2]=eig(K2);
    d1=sort(diag(real(D1)),'descend'); d2=sort(diag(real(D2)),'descend');
    %figure(); plot(d1); hold on; plot(d2);
    Dim=min(find(d1<th,1,'first'),find(d2<th,1,'first'));
    if isempty(Dim), Dim=-1;,end;
end

