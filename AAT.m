function y = AAT(x,flag)
global A;
global N_dim;

if nargin <2
y = A*(A'*x);

elseif strcmp(flag,'dimension')
    y = N_dim;
    
else
   y = [];
end

return