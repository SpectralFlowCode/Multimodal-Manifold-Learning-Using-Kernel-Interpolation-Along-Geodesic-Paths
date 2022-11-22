N=DataParams.NumberOfDatapoints;

if strcmp(DataParams.Layout,'Torus_Torodial') %common is dominant
    X=2*pi*(rand(1,N)-0.5);
    Y=2*pi*(rand(1,N)-0.5);
    Z=2*pi*(rand(1,N)-0.5);
    
    X(1)=pi+pi/2;Y(1)=pi;Z(1)=pi;
    S1=[(DataParams.R1+DataParams.r1*cos(Y)).*cos(X);
        (DataParams.R1+DataParams.r1*cos(Y)).*sin(X);
        DataParams.r1*sin(Y)]';
    S2=[(DataParams.R2+DataParams.r2*cos(Z)).*cos(X);
        (DataParams.R2+DataParams.r2*cos(Z)).*sin(X);
        DataParams.r2*sin(Z)]';
    
    
end
if strcmp(DataParams.Layout,'Torus_Polodial') %common is small (common is the polodial angle)
    X=2*pi*(rand(1,N)-0.5);
    Y=2*pi*(rand(1,N)-0.5);
    Z=2*pi*(rand(1,N)-0.5);
    
    X(1)=0;Y(1)=pi;Z(1)=pi;
    S1=[(DataParams.R1+DataParams.r1*cos(X)).*cos(Y);
        (DataParams.R1+DataParams.r1*cos(X)).*sin(Y);
        DataParams.r1*sin(X)]';
    S2=[(DataParams.R2+DataParams.r2*cos(X)).*cos(Z);
        (DataParams.R2+DataParams.r2*cos(X)).*sin(Z);
        DataParams.r2*sin(X)]';
end
if strcmp(DataParams.Layout(1:6),'2DFlat')   %2D flat manifold 
    X=2*(rand(1,N)-0.5);
    Y=2*(rand(1,N)-0.5);
    Z=2*(rand(1,N)-0.5);
    
    X(1)=0;Y(1)=0;Z(1)=0;
    S1=[DataParams.ScaleX*X;DataParams.ScaleY*Y]';
    S2=[DataParams.ScaleX*X;DataParams.ScaleZ*Z]';    
end


