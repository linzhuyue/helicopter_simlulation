function flag = ssorder(s,v,tol)

%SSORDER  Ordering of Vector Subspace
%
%         ss = ssorder(X,Y)
%
%         determines the ordering of two vector spaces respectively
%         spanned by the columns of matrices X and Y.
%
%         Output Parameters:
%
%          if ss = -1, subspace spanned by X < that spanned by Y
%          if ss =  0, subspace spanned by X = that spanned by Y
%          if ss =  1, subspace spanned by X > that spanned by Y
%          if ss =  j = sqrt(-1), they are not related at all.

if nargin==2
   tol=1e-8;
end;

ns=size(s,1);
nv=size(v,1);
%if (ns~=nv)&(ns~=0)&(nv~=0)
if (ns~=nv)
   flag=j;return;
end

ns=size(zzorthtol(s,tol),2);
nv=size(zzorthtol(v,tol),2);
nsv=size(zzorthtol([s,v],tol),2);

if nsv>max(ns,nv)
   flag=j;
end
if nsv==ns
   flag=1;
end
if nsv==nv
   flag=-1;
end
if (nsv==ns) & (nsv==nv)
   flag=0;
end