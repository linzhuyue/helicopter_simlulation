function S=r_star(A,B,C,D,tol)

%R_STAR  Controllable Weakly Unobservable Geometric Subspace
%
%        R = r_star(A,B,C,D)
%
%        computes a matrix whose columns span the geometric subspace
%        R^{*} for a system characterized by (A,B,C,D).
%
%        Note that this function is applicable for both continuous-
%        and discrete-time systems.
%
%        See also S_STAR, V_STAR and N_STAR.

if nargin==4,
   tol=1e-8;
end

flag='r_star';
dc=0;
S=zzscb_space(A,B,C,D,flag,dc,tol);

if isempty(S)
   S=zeros(size(A,1),1);
else
   S=orth(S);
end