function VV=ssintsec(s,v,tol)

%SSINTSEC  Intersection of Vector Subspace
%
%          V = ssintsec(X,Y)
%
%          computes intersection of two vector spaces respectively
%          spanned by the columns of matrices X and Y.
%
%          The columns of V form a basis for the intersection.
%
%          See also SSORDER and SSADD.


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

if isempty(s)|isempty(v)
   VV=[];
   return
end

if nargin==2
   s=null(s');
   v=null(v');
   VV=ssadd(s,v);
   VV=null(VV');
else
   s=zznulltol(s',tol);
   v=zznulltol(v',tol);   
   VV=ssadd(s,v,tol);
   VV=zznulltol(VV',tol);
end

if isempty(VV)
   VV=zeros(nv,1);
end