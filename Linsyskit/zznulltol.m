function Z = zznulltol(A,tol)
%NULL   Null space.
%   Z = zznulltol(A) is an orthonormal basis for the null space of A obtained
%   from the singular value decomposition.  That is,  A*Z has negligible
%   elements, size(Z,2) is the nullity of A, and Z'*Z = I.
%

if nargin==1
   tol=1e-10;
end

[m,n] = size(A);
[U,S,V] = svd(A,0);
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
Z = V(:,r+1:n);