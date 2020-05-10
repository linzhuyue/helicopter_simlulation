function I2=r_invt(A,B,C,D,tol);

%R_INVT  Right Invertibility Structure of Proper Systems
%
%        rights = r_invt(A,B,C,D)
%
%        returns the right invertibility structure of a multivariable
%        system characterized by (A,B,C,D).
%
%           rights = I_2 List of Morse Indices
%
%        See also INVZ, INFZ, R_INVT and MORSEIDX.

if nargin==4,
   tol=1e-8;
end

dc=0;
[AA,BB,CC,DD,Gs,Go,Gi,dims,t,I2]=scb(A,B,C,D,tol,dc);

%type=2; dc=0; d11_eye=0;
%[AA,BB,CC,DD,Gs,Go,Gi,dims,I3,I2,I4,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);
