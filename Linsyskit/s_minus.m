function S=s_minus(A,B,C,D,dc,tol)

%S_MINUS  Stable Strongly Controllable Geometric Subspace
%
%         S = s_minus(A,B,C,D[,dc])
%
%         computes a matrix whose columns span the geometric subspace
%         S^{-} for a system characterized by (A,B,C,D).
%
%         Note that by default or if dc = 0, the function returns a
%         subspace for a continuous-time system. Otherwise, if
%         dc = 1, it computes a subspace for a discrete-time system.
%
%         See also S_STAR, S_PLUS and V_MINUS.

if nargin==4,
   dc=0; tol=1e-8;
elseif nargin==5
   tol=1e-8;
end

flag='s_minu';
S=zzscb_space(A,B,C,D,flag,dc,tol);

if isempty(S)
   S=zeros(size(A,1),1);
else
   S=orth(S);
end