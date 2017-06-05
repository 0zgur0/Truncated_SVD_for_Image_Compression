function [U,Sigma,V] = trunc_svd2b(A,r) 
%Method 2, large-scale SVD

[m,n] = size(A);

%Check if r is valid 
if r>=min(m,n)
    display('r has to be less than or eqaul to min(m,n), r is set to min(m,n)');
    r = min(m,n);
end

%Define alpha
alpha = norm(A)/10;

%if m>=n, use A'A otherwise AA'
if m>=n

    M1 = chol(A'*A+alpha*eye(n),'lower'); M2 = M1';
    %Precond matrix
    M = [M1,M2];

    [V,D] = jdqr(A'*A,'K',r,struct('Precond',M));
    Sigma = sqrt(D);
    temp = diag(sqrt(D)); 
    Sigma_inv = diag(1./temp);
    U = A*V*Sigma_inv;  

    
else
    
    M1 = chol(A*A'+alpha*eye(m),'lower'); M2 = M1';
    %Precond matrix
    M = [M1,M2];
    
    [U,D] = jdqr(A*A','K',r,struct('Precond',M));
    Sigma = sqrt(D);
    
    temp = diag(sqrt(D));
    Sigma_inv = diag(1./temp);
    V = Sigma_inv*(U')*A;
    V = V';
    
    
end



end