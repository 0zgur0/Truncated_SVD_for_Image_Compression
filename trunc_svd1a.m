function [U,Sigma,V] = trunc_svd1a(A,r) 
%Method 1a, thin SVD
%This function computes the truncated SVD of matrix A
%r is desired number of the most dominant eigenvalues
% r has to be less than or eqaul to min(m,n)
% U is m by r with orthonormal columns
%V is n by r with orthonormal columns
%Sigma is r by r diagonal matrix with dominant eigenvalues

[m,n] = size(A);

%Check if r is valid 
if r>=min(m,n)
    display('r has to be less than or eqaul to min(m,n), r is set to min(m,n)');
    r = min(m,n);
end

[U,Sigma,V] = svd(A,'econ');
U = U(:,1:r);
V = V(:,1:r);
Sigma = Sigma(1:r,1:r);

end

