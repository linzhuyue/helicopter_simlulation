function [Zeig,I1]=invz(A,B,C,D,tol);

%INVZ  Invariant Zeros and Structures of Proper Systems
%
%      zrs = invz(A,B,C,D)
%
%      returns invariant zeros of a system characterized by (A,B,C,D)
%      and their structures.
%
%         zrs = all the invariant zeros of the system
%
%      Note that invariant zeros are sometimes called transmission
%      zeros. However, the latter is only defined for controllable
%      and observable systems.
%
%      See also BLKZ and INFZ and MORSEIDX.

if nargin==4,
   tol=1e-8;
end

dc=0;
[AA,BB,CC,DD,Gs,Go,Gi,dims]=scb(A,B,C,D,tol,dc);

Zeig=[];
at1=AA(1:dims(1),1:dims(1));
[I1,T,struc,Zeig]=jcf(at1,tol);
if ~isempty(Zeig)
   Zeig=Zeig(:,1)+i*Zeig(:,2);
end

at1=AA(dims(1)+1:sum(dims(1:2)),dims(1)+1:sum(dims(1:2)));
[jt,T,struc,Zt]=jcf(at1,tol);
I1=blkdiag(I1,jt);
if ~isempty(Zt)
   Zeig=[Zeig;Zt(:,1)+i*Zt(:,2)];
end

at1=AA(sum(dims(1:2))+1:sum(dims(1:3)),sum(dims(1:2))+1:sum(dims(1:3)));
[jt,T,struc,Zt]=jcf(at1,tol);
I1=blkdiag(I1,jt);
if ~isempty(Zt)
   Zeig=[Zeig;Zt(:,1)+i*Zt(:,2)];
end

Zeig=sort(eig(AA(1:sum(dims(1:3)),1:sum(dims(1:3)))));