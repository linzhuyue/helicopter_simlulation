function Q = zzorthtol(A,tol)
%ORTH   Orthogonalization.
%   Q = zzorthtol(A) is an orthonormal basis for the range of A.
%   That is, Q'*Q = I, the columns of Q span the same space as 
%   the columns of A, and the number of columns of Q is the 
%   rank of A.
%

if nargin==1
   tol=1e-10;
end

[U,S,V] = svd(A,0);
[m,n] = size(A);
if m > 1
   s = diag(S);
%elseif m == 1
elseif (m==1)&(n~=0)
   s = S(1);
else
   s = 0;
end
%tol = max(m,n) * max(s) * eps;
r = sum(s > tol);
Q = U(:,1:r);