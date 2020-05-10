function I4=infz(A,B,C,D,tol);

%INFZ  Infinite Zero Structure of Proper Systems
%
%      infzs = infz(A,B,C,D)
%
%      returns the infinite zero structure of a system characterized
%      by (A,B,C,D).
%
%         infzs = I_4 List of Morse Indices
%
%      See also INVZ, L_INVT, R_INVT and MORSEIDX.

if nargin==4,
   tol=1e-8;
end

dc=0;
[AA,BB,CC,DD,Gs,Go,Gi,dims,t,t,I4]=scb(A,B,C,D,tol,dc);

%type=1; dc=0; d11_eye=0;
%[AA,BB,CC,DD,Gs,Go,Gi,dims,I3,I2,I4,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);
