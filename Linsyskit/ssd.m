function      [AA,T,nn,no,np,err_SSD]=ssd(A,tol)

%SSD  Stability Structural Decomposition 
%
%     [D,T,nn,no,np] = ssd(A[,tol])
%
%     gives the following block diagonal form for a square matrix:
%
%                         [ A-   0    0  ]  nn
%        D = inv(T)*A*T = | 0    Ao   0  |  no
%                         [ 0    0    A+ ]  np
%
%     where eigenvalues of A-, Ao and A+ are, respectively, in the
%     open left-half plane, jw axis and open right-half plane.
%
%     See also DSSD.

if nargin==1
   tol=1e-8;
end;

flag=1;
[AA, T, nn, no, np, err_SSD]=zzcssdresch(A,tol);
[AA2,T2,nn2,no2,np2,err2]=zzcssdblkr(A,tol);

if err2<err_SSD
   AA=AA2;
   T=T2;
   nn=nn2;
   no=no2;
   np=np2;
   err_SSD=err2;
end
