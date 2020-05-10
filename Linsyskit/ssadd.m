function VV=ssadd(s,v,tol)

%SSADD  Addition of Vector Subspace
%
%       V = ssadd(X,Y)
%
%       computes the addition of two vector spaces respectively
%       spanned by the columns of matrices X and Y.
%
%       The columns of V form a basis for the addition.
%
%       See also SSORDER and SSINTSEC.

ns=size(s,1);
nv=size(v,1);
%if (ns~=nv)&(ns~=0)&(nv~=0)
if (ns~=nv)
   disp(' ')
   disp(' The rows of X and Y are different, ... ...')
   disp(' ')
   VV=[];
   return
end

if isempty(s)&isempty(v)
   VV=s;
   return;
end

if nargin==2
   VV=orth([s,v]);
else
   VV=zzorthtol([s,v],tol);
end

if isempty(VV)
   VV=zeros(nv,1);
end