function [AA,T,struc,err_zzbdeig]=zzbdeig(A,tol)

%ZZBDEIG
%
%     [Aj,P,zstruc]=zzbdeig(A[,tol])
%
%     transforms a square matrix A into a block-diagonalize matrix whose blocks
%     have either a single or one repeated real eigenvalue, or a single or one
%     repeated pair of complex eigenvalues.
%
%   See also CSSD.

if nargin==1,
   tol=1e-6;
end;

n=size(A,1);T=eye(n);
Akk=A;
struc=[];
while ~isempty(Akk)
   eA=sort(real(eig(Akk,'nobalance')));
   mu1=min(eA);
   n1kk=size(Akk,1);
   kk=sum(eA<(mu1+2*tol));
   if kk<n1kk
      muk1=eA(kk+1);
      Abar=Akk-(mu1+muk1)/2*eye(n1kk);
      [At,T1,nn,no,np]=ssd(Abar,tol);
      At=inv(T1)*Akk*T1;
      tt=eye(n);tt(n-n1kk+1:n,n-n1kk+1:n)=T1;
      T=T*tt;
      A1k=At(1:kk,1:kk);
      Akk=At(kk+1:n1kk,kk+1:n1kk);
   else
      A1k=Akk;
      Akk=[];
   end
   
   nn1k=size(A1k,1);T1k=eye(nn1k);
   while ~isempty(A1k)
      n1k=size(A1k,1);
      eAk=sort(abs(imag(eig(A1k,'nobalance'))));
      omiga1=min(eAk);
      k1=sum(eAk<(omiga1+2*tol));
      struc=[struc,k1];
      if k1<n1k
         omiga2=eAk(k1+1);
         beta=sqrt(mu1^2+omiga1^2)/2+sqrt(mu1^2+omiga2^2)/2;
         Ahat=inv(A1k+beta*eye(n1k))*(A1k-beta*eye(n1k));
         [At,T2,nn,no,np]=ssd(Ahat,tol);
         At=inv(T2)*A1k*T2;
         tt=eye(nn1k);tt(nn1k-n1k+1:nn1k,nn1k-n1k+1:nn1k)=T2;
         T1k=T1k*tt;
         A1k=At(k1+1:n1k,k1+1:n1k);
      else
         A1k=[];
      end
   end
   
   ttt=eye(n1kk);ttt(1:nn1k,1:nn1k)=T1k;
   tt=eye(n);tt(n-n1kk+1:n,n-n1kk+1:n)=ttt;
   T=T*tt;
end

if cond(T)>1e15
   AA=A;T=eye(n);struc=[];err_zzbdeig=inf;flag=0;
   return
end

% Makes the elments should be zero be zero.
tt=[];
for k=1:length(struc)
    tt=blkdiag(tt,ones(struc(k)));
end
AA=inv(T)*A*T.*tt;
err_zzbdeig=norm(AA-inv(T)*A*T);