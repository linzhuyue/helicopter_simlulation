function qv = infz_ds(E,A,B,C,D,tol)

%INFZ_DS  Infinite Zero Structure of Descriptor Systems
%
%         infzs = infz_ds(E,A,B,C,D)
%
%         returns the infinite zero structure of a descriptor
%         system characterized by (E,A,B,C,D).
%
%         See also INFZ.

if nargin==5
   tol=1e-5;
end

[Es4,As4,Bs4,Cs4,Ds4,Ez4,B0,C0,L4,Ck4,Dk4,yrb4,Ge,Gs,Go,invGi,dims,lv,rv,qv,m0]=sd_ds(E,A,B,C,D,tol);
