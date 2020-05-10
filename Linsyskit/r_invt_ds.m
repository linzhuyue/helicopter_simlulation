function rv = r_invt_ds(E,A,B,C,D,tol)

%R_INVT_DS  Right Invertibility Structure of Descriptor Systems
%
%           rights = r_invt_ds(E,A,B,C,D)
%
%           gives the right invertibility structure of a descriptor
%           system characterized by (E,A,B,C,D).
%
%           See also R_INVT.

if nargin==5
   tol=1e-5;
end

[Es4,As4,Bs4,Cs4,Ds4,Ez4,B0,C0,L4,Ck4,Dk4,yrb4,Ge,Gs,Go,invGi,dims,lv,rv,qv,m0]=sd_ds(E,A,B,C,D,tol);