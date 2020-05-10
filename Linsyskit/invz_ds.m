function [Zeig,I1] = invz_ds(E,A,B,C,D,tol)

%INVZ_DS  Invariant Zeros and Structures of Descriptor Systems
%
%         zrs = invz_ds(E,A,B,C,D)
%
%         gives invariant zeros of a descriptor system characterized
%         by (E,A,B,C,D) and their structures.
%
%           zrs = all the invariant zeros of the system
%
%         See also INVZ.

if nargin==5
   tol=1e-5;
end

[Es4,As4,Bs4,Cs4,Ds4,Ez4,B0,C0,L4,Ck4,Dk4,yrb4,Ge,Gs,Go,invGi,dims,lv,rv,qv,na,m0]=sd_ds(E,A,B,C,D,tol);

AA=As4(sum(dims(1:2))+1:sum(dims),sum(dims(1:2))+1:sum(dims));

dims=na;

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