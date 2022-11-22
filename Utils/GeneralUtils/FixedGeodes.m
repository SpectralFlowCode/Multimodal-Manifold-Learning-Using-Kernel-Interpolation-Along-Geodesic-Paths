function [ S ] = FixedGeodes( A,B,t,Dim )
% Input:
%   A,B - PSD matrices of rank Dim
%   0<t<1 - Desired point along the geodesic
%   Dim - rank of A or B.
% Output:
%   S - the point t along the geodesic streching from A to B.

if Dim>0 
    % Apply fixed rank approximation according to: Bonnabel, Silvere, and
    % Rodolphe Sepulchre. "Riemannian metric and geometric mean for positive semidefinite matrices of fixed rank."
    % SIAM Journal on Matrix Analysis and Applications 31.3 (2010): 1055-1070.
    
    [U1,~]=eigs(A,Dim);
    [U2,~]=eigs(B,Dim);
      
    VA=U1(:,1:Dim);VB=U2(:,1:Dim);
    [OA,SAB,OB]=svd(VA'*VB);
    % Allgedly the following is redundant, however, it seems to be more
    % stable as it mangae to fix numerical errors:
    SAB(SAB>1)=1;
    UA=VA*OA;UB=VB*OB;
    theta=acos(diag(SAB));Theta=diag(theta);
    X=(eye(size(A))-UA*transpose(UA))*UB*pinv(diag(sin(theta)));
    U=UA*diag(cos(theta*t))+X*diag(sin(theta*t));
    
    RA2=transpose(UA)*A*UA;%RA2=0.5*(RA2+RA2');
    RB2=transpose(UB)*B*UB;%RB2=0.5*(RB2+RB2');
    %     Suppose to be:
    %     RA=(RA2)^(1/2);RB=(RB2)^(1/2);
    %     Instead, I use the following computation, which is more stable:
    [UR1,SR1,VR1]=svd(RA2,0);
    RA=UR1*(SR1.^(1/2))*VR1';
    [UR2,SR2,VR2]=svd(RB2,0);
    RB=UR2*(SR2.^(1/2))*VR2';
    R2=RA*expm(t*logm(RA^(-1)*RB^2*RA^(-1)))*RA;
    S=U*R2*transpose(U);
    
else
    %Apply full rank computation
    switch t
        case 0
            S=A;
        case 1
            S=B;
        otherwise
            % Instead of using:
            %     S=(A^(0.5))*(((A^(-0.5))*B*(A^(-0.5)))^t)*(A^(0.5));
            % We propose using the following, which is more numerically stable:
            [U,S,V] = svd(abs((A^(-0.5))*B*(A^(-0.5))));
            S=(A^(0.5))*(U*S.^t*V')*(A^(0.5));
    end
    
end


end



