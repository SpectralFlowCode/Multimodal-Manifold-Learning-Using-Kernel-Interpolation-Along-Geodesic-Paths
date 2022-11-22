function  [R2,yest]= GetPolyFitLoss(x,y,p)
    pol = polyfit(x,y,p);
    yest = polyval(pol,x);
    R2=norm(yest-y)/norm(y);
end

