function NR = normrank(A,B,C,D,tol)

%NORMRANK  Normal Rank of Proper System
%
%          NR = normrank(A,B,C,D)
%
%          returns the normal rank of a linear system characterized
%          by (A,B,C,D).
%
%          See also INVZ.

if nargin==4,
   tol=1e-8;
end

dc=0;
[AA,BB,CC,DD,Gs,Go,Gi,dims,t,t,I4,m0]=scb(A,B,C,D,tol,dc);

%type=1; dc=0; d11_eye=0;
%[AA,BB,CC,DD,Gs,Go,Gi,dims,I3,I2,I4,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);

NR=length(I4)+m0;