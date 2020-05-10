function S=v_plus(A,B,C,D,dc,tol)

%V_PLUS  Unstable Weakly Unobservable Geometric Subspace
%
%        V = v_plus(A,B,C,D[,dc])
%
%        computes a matrix whose columns span the geometric subspace
%        V^{+} for a system characterized by (A,B,C,D).
%
%        Note:
%
%        By default or if dc = 0, the function returns a subspace for
%        a continuous-time system. Otherwise, if dc = 1, it computes
%        a subspace for a discrete-time system.
%
%        See also V_STAR, V_MINUS and S_PLUS.

if nargin==4,
   dc=0; tol=1e-8;
elseif nargin==5
   tol=1e-8;
end

flag='v_plus';
S=zzscb_space(A,B,C,D,flag,dc,tol);

if isempty(S)
   S=zeros(size(A,1),1);
else
   S=orth(S);
end