function [X,flag] = zzare(A,B,C)

%ARE  Algebraic Riccati Equation solution.
%
%    Note that the function is rewritten from are.m in Control toolbox.

%  -- check for correct input problem --
[nr,nc] = size(A); n = nr;
if (nr ~= nc), error('Nonsquare A matrix'); end;
[nr,nc] = size(B);
if (nr~=n | nc~=n), error('Incorrectly dimensioned B matrix'); end;
[nr,nc] = size(C);
if (nr~=n | nc~=n), error('Incorrectly dimensioned C matrix'); end;

% Following is much faster than before
[q,t] = schur([A -B; -C -A']);
[q,t] = rsf2csf(q,t);

tol = 10.0*eps*max(abs(diag(t)));   % ad hoc tolerance
ns = 0;
%
%  Prepare an array called index to send message to ordering routine 
%  giving location of eigenvalues with respect to the imaginary axis.
%  -1  denotes open left-half-plane
%   1  denotes open right-half-plane
%   0  denotes within tol of imaginary axis
%
index = [];  
for i = 1:2*n,
    if (real(t(i,i)) < -tol),
        index = [ index -1 ];
    ns = ns + 1;
    elseif (real(t(i,i)) > tol),
    index = [ index 1 ];
    else,
    index = [ index 0 ];
    end;
end;
flag=1;
if (ns ~= n),
   flag=0;
   X=[];
   return
end
[q,t] = schord(q,t,index);
X = real(q(n+1:n+n,1:n)/q(1:n,1:n));

if any(any(X==inf))
   flag=0;
else
   flag=(min(real(eig(X)))>tol);
end