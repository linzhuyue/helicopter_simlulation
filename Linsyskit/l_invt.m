function I3=l_invt(A,B,C,D,tol);

%L_INVT  Left Invertibility Structure of Proper Systems
%
%        lefts = l_invt(A,B,C,D)
%
%        returns the left invertibility structure of a multivariable
%        system characterized by (A,B,C,D).
%
%           lefts = I_3 List of Morse Indices
%
%        See also INVZ, INFZ, R_INVT and MORSEIDX.

if nargin==4,
   tol=1e-8;
end

dc=0;
[AA,BB,CC,DD,Gs,Go,Gi,dims,I3]=scb(A,B,C,D,tol,dc);

%type=2; dc=0; d11_eye=1;
%[AA,BB,CC,DD,Gs,Go,Gi,dims,I3,I2,I4,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);
