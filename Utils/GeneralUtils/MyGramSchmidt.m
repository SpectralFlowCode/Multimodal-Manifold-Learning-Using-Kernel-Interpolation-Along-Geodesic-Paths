function OV = MyGramSchmidt(V)

OV=[V(:,1)/norm(V(:,1))];
for ind=2:size(V,2)
    tmpV=V(:,ind);%/norm(V(:,ind));
    set_projection=tmpV*0;
    for sub_ind=1:size(OV,2)
        set_projection=set_projection+OV(:,sub_ind)'*tmpV*OV(:,sub_ind);
    end
    tmpV=tmpV-set_projection;
    tmpV=tmpV/norm(tmpV);
    OV=[OV,tmpV];
end

end