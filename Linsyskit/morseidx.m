function [I1,I2,I3,I4]=morseidx(A,B,C,D,tol);

%MORSEIDX  Morse Invariance Indices of Proper Systems
%
%          [I1,I2,I3,I4] = morseidx(A,B,C,D)
%
%          returns Morse structural invariance list sfor a system
%          characterized by (A,B,C,D).
%
%            I1 = zero dynamics matrix in Jordan form
%            I2 = I_2 List = right invertibility structure
%            I3 = I_3 List = left invertibility structure
%            I4 = I_4 List = infinite zero structure
%
%          Note that I1 List should formally contain the invariant
%          zeros and the sizes of their Jordan blocks.
%
%          See also INVZ, INFZ, L_INVT and R_INVT.

if nargin==4,
   tol=1e-8;
end

dc=0;
[AA,BB,CC,DD,Gs,Go,Gi,dims,I3,I2,I4]=scb(A,B,C,D,tol,dc);

%type=2; dc=0; d11_eye=0;
%[AA,BB,CC,DD,Gs,Go,Gi,dims,I3,I2,I4,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);

at1=AA(1:dims(1),1:dims(1));
I1=jcf(at1,tol);
at1=AA(dims(1)+1:sum(dims(1:2)),dims(1)+1:sum(dims(1:2)));
jt=jcf(at1,tol);
I1=blkdiag(I1,jt);
at1=AA(sum(dims(1:2))+1:sum(dims(1:3)),sum(dims(1:2))+1:sum(dims(1:3)));
jt=jcf(at1,tol);
I1=blkdiag(I1,jt);
