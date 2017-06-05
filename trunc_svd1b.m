function [U,Sigma,V] = trunc_svd1b(A,r)
%Method 1b, thin QR, then SVD

[m,n] = size(A);

%Check if r is valid
if r>=min(m,n)
    display('r has to be less than or eqaul to min(m,n), r is set to min(m,n)');
    r = min(m,n);
end

if m>=n
    %QR factorization of input matrix A
    %Option '0' enables thin QR factorization
    [Q,R] = qr(A,0);
    
    %SVD of square matrix R
    [U,Sigma,V] = svd(R);
    
    U = Q*U;
    
    
else % in the case of m<n, then we should work with A^T
    %QR factorization of transpose of input matrix A
    [Q,R] = qr(A',0);
    
    %SVD of square matrix R
    [U1,Sigma,V1] = svd(R);
    
    
    U = V1;
    V = Q*U1;
    
end

U = U(:,1:r);
V = V(:,1:r);
Sigma = Sigma(1:r,1:r);

end