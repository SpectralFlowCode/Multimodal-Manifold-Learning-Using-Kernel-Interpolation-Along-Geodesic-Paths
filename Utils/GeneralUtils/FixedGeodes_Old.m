function [ S ] = FixedGeodes( A,B,t,Dim )
% Input:
%   A,B - PSD matrices of rank Dim
%   0<t<1 - Desired point along the geodesic
%   Dim - rank of A or B.
% Output:
%   S - the point t along the geodesic streching from A to B.
if Dim>0
    
    [U1,S1]=eigs(A,Dim);%
    [U2,S2]=eigs(B,Dim);%
    %figure; ax(1)=subplot(2,1,1);plot(sum(abs(U1-V1))); ax(2)=subplot(2,1,2);plot(db(diag(S1))); linkaxes(ax,'x');
    
    
    VA=U1(:,1:Dim);VB=U2(:,1:Dim);
    [OA,SAB,OB]=svd(VA'*VB);
%     SAB(SAB<1)=1;
    UA=VA*OA;UB=VB*OB;
    theta=acos(diag(SAB));Theta=diag(theta);
    X=(eye(size(A))-UA*transpose(UA))*UB*pinv(diag(sin(theta)));
    U=UA*diag(cos(theta*t))+X*diag(sin(theta*t));
    
    RA2=transpose(UA)*A*UA;%RA2=0.5*(RA2+RA2');
    RB2=transpose(UB)*B*UB;%RB2=0.5*(RB2+RB2');
    %     Suppose to be:
    RA=(RA2)^(1/2);RB=(RB2)^(1/2);
    %     Instead, I use:
    %     RA=sqrt(S1(1:Dim,1:Dim));RB=sqrt(S2(1:Dim,1:Dim));
    R2=RA*expm(t*logm(RA^(-1)*RB^2*RA^(-1)))*RA;
    S=U*R2*transpose(U);
    
else
    switch t
        case 0
            S=A;
        case 1
            S=B;
        otherwise
            % Instead of using:
            %     S=(A^(0.5))*(((A^(-0.5))*B*(A^(-0.5)))^t)*(A^(0.5));
            % we propose using the following, which is more numerically stable:
            tmp=abs((A^(-0.5))*B*(A^(-0.5)));
            %tmp=abs((pinv(A)^(0.5))*B*(pinv(A)^(0.5)));
            [U,S,V] = svd(tmp);
            tmp2= U*S.^t*V';
            S=(A^(0.5))*tmp2*(A^(0.5));
    end
    
end


end



