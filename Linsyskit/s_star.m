function S=s_star(A,B,C,D,tol)

%S_STAR  Strongly Controllable Geometric Subspace
%
%        S = s_star(A,B,C,D)
%
%        computes a matrix whose columns span the geometric subspace
%        S^{*} for a system characterized by (A,B,C,D).
%
%        Note that this function is applicable for both continuous-
%        and discrete-time systems.
%
%        See also S_MINUS, S_PLUS and V_STAR.

if nargin==4,
   tol=1e-8;
end;

flag='s_star';
dc=0;
S=zzscb_space(A,B,C,D,flag,dc,tol);

if isempty(S)
   S=zeros(size(A,1),1);
else
   S=orth(S);
end