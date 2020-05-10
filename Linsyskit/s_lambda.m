function S=s_lambda(A,B,C,D,lambda,tol)

%S_LAMBDA  Geometric Subspace S_lambda (See Chapter 3 for definition)
%
%          S = s_lambda(A,B,C,D,lambda)
%
%          computes a matrix whose columns span S_{\lambda} for a
%          system characterized by (A,B,C,D), where 'lambda' is
%          either a real or complex scalar.
%
%          Note: This function is applicable for both continuous- and
%          discrete-time systems.
%
%          See also V_LAMBDA.

if nargin==5
   tol=1e-8;
end;

flag='s_lamb';
dc=0;
S=zzscb_space(A,B,C,D,flag,dc,tol,lambda);

if isempty(S)
   S=zeros(size(A,1),1);
else
   S=orth(S);
end