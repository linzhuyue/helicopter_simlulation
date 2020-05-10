function [AA,T,nn,no,np,err_of_DSSD]=zzdssd3s(A,tol)

%DSSD3S  Stability Structural Decomposition for Discrete-time Systems
%
%     [A,T,nn,no,np]=zzdssd3s(A[,tol])
%
%     puts any square matrix in the following block diagonal form,
%
%                               | A-    0     0  |  nn
%             AA = inv(T)*A*T = | 0     Ao    0  |  no
%                               | 0     0     A+ |  np
%
%     where eigenvalues of (A-), (Ao) and (A+) are, respectively,
%     inside, on and outside the unit circle of the complex plane.
%
%     See also DSSD, DSSDalpha, SSD.

%    Note that: DSSD3S uses the function cssd for three times.

%     Ben M. Chen
%     National University of Singapore
%     March 27, 1994; Revised October 5, 1996

if nargin==1
   tol=1e-10;
end;

n=size(A,1);

[AA,T,nn,no,np]=ssd(A,tol);

AA=inv(T)*A*T;

Tno=[];Tp=[];

if n==0
   AA=[];T=[];nn=0;no=0;np=0;err_of_DSSD=0;
   return;
end

if nn+no>0
   Ax=AA(1:nn+no,1:nn+no);
   I=eye(nn+no);
   At=inv(Ax-I)*(Ax+I);
   [AAno,Tno,nnno,nono,npno]=ssd(At,tol);
else
   nnno=0;nono=0;npno=0;
end

if np>0
   Ax=AA(nn+no+1:n,nn+no+1:n);
   I=eye(np);
   At=inv(Ax+I)*(Ax-I);
   [AAp,Tp,nnp,nop,npp]=ssd(At,tol);
else 
   nnp=0;nop=0;npp=0;
end

t=eye(n);
Tq=[t(1:nnno,:);t(nn+no+1:nn+no+nnp,:);t(nnno+1:nnno+nono,:); ......
      t(nn+no+nnp+1:nn+no+nnp+nop,:);t(nnno+nono+1:nnno+nono+npno,:);t(n-npp+1:n,:)];
T=T*[Tno zeros(nn+no,np);zeros(np,nn+no) Tp]*inv(Tq);
nn=nnno+nnp;
no=nono+nop;
np=npno+npp;

AA=inv(T)*A*T;
tt=blkdiag(ones(nn,nn),ones(no,no),ones(np,np));
AA=AA.*tt;

err_of_DSSD=norm(AA-inv(T)*A*T);