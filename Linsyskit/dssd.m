function [AA,T,nn,no,np,err_of_DSSD]=dssd(A,tol)

%DSSD  Discrete Version Stability Structural Decomposition
%
%      [D,T,nn,no,np] = dssd(A[,tol])
%
%      gives the following block diagonal form for a square matrix:
%
%                          [ A-   0    0  ]  nn
%         D = inv(T)*A*T = | 0    Ao   0  |  no
%                          [ 0    0    A+ ]  np
%
%      where eigenvalues of A-, Ao and A+ are, respectively, inside,
%      on and outside the unit circle of the complex plane.
%
%      See also SSD.
 
if nargin==1
   tol=1e-8;
end

[AA, T, nn, no, np, err_of_DSSD]=zzdssd3s(A,tol);
[AA2,T2,nn2,no2,np2,errm]=zzdssdalpha(A,tol);

if errm<err_of_DSSD
   AA=AA2;
   T=T2;
   nn=nn2;
   no=no2;
   np=np2;
   err_of_DSSD=errm;
end