function      [AA,T,nn,no,np,errm]=zzdssdalpha(A,tol)

%DSSDALPHA  Stability Structural Decomposition for Discrete-time Systems
%
%     [A,T,nn,no,np]=zzdssdalpha(A[,tol])
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
%    See also DSSD, DSSD3S, SSD.

if nargin==1
   tol=1e-12;
end;

n=size(A,1);
if n==0
   AA=[];T=[];nn=0;no=0;np=0;errm=0;
   return;
end

eA=eig(A);errm=inf;
for kk=1:3
   mt=0;
   for k=1:5
      b1=2*pi*rand; b1=sin(b1)+i*cos(b1);
      tt=[abs(eA-b1),abs(eA+b1),abs(eA-conj(b1)),abs(eA+conj(b1))];
      tt=sort(tt(:)); t1=tt(1); t2=sum(tt(1:4)); t3=1*t1+1*t2;
      if t3>mt
         mt=t3;
         alpha=b1;
      end
   end
   Abar=(inv(A+alpha*eye(n))*(A-alpha*eye(n))+inv(A+conj(alpha)*eye(n))*(A-conj(alpha)*eye(n)))/2;
   Abar=real(Abar);
   [AA,T,nn,no,np]=ssd(Abar,tol);
   AA=inv(T)*A*T;
   tt=blkdiag(ones(nn,nn),ones(no,no),ones(np,np));
   AA=AA.*tt;   
   et=norm(AA-inv(T)*A*T);
   if et<errm
      errm=et;
      Tm=T;Am=AA;
   end   
end
T=Tm;AA=Am;
%errm