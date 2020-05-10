function S=v_star(A,B,C,D,tol)

%V_STAR  Weakly Unobservable Geometric Subspace
%
%        V = v_star(A,B,C,D)
%
%        computes a matrix whose columns span the geometric subspace
%        V^{*} for a system characterized by (A,B,C,D).
%
%        Note:
%
%        It is applicable for both continuous- and discrete-time
%        systems.
%
%        See also V_MINUS, V_PLUS and S_STAR.

if nargin==4,
   tol=1e-8;
end

flag='v_star';
dc=0;
S=zzscb_space(A,B,C,D,flag,dc,tol);

if isempty(S)
   S=zeros(size(A,1),1);
else
   S=orth(S);
end