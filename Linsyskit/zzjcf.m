function [AA,T,struc,Jeig,err_of_JCF,flag]=zzjcf(A,tol)

%ZZJCF  Jordan Canonical Form
%
%     [J,T] = zzjcf(A[,tol])
%
%     generates a transformation that transforms a square matrix into
%     the usual Jordan canonical form, i.e.,
%
%                           [ J_1         ]
%         inv(T)*A*T = J =  |     .       |
%                           |       .     |
%                           [         J_k ]
%
%      where each block J_i, i=1, ..., k, has the following form:
%
%                      [ eig_i    1          ]
%                      |     .       .       |
%               J_i =  |       .       .     |
%                      |         eig_i   1   |
%                      [               eig_i ]
%
%      See also RJD.

%     stru(i,:)= the sizes of jordan blocks with the eigenvalue mu_i + j * omega_i;
%     Jeig(i,:)=[mu_i, omega_i];

%     flag=1, the decomposition is reasonable for the given tol.
%     flag=0, the decomposition can not be carried out for the given tol.


if nargin==1
   tol=1e-6;
end

n=size(A,1);
if n==0
   AA=[];T=[];struc=0;Jeig=[];err_of_JCF=[];flag=1;
   return;
end

[AA,T,st,RJ,err_of_RJD,tolt,flag]=rjd(A,tol);
if flag==0
   AA=A;T=eye(n);struc=[];Jeig=[];err_of_JCF=inf;
   return;
end

flag=1;
struc=[];Jeig=zeros(0,2);TT=[];
for kk=1:size(st,1)
   t1=sum(sum(st(1:kk-1,:))); t2=sum(st(kk,:));
   if RJ(kk,2)==0
      TT=blkdiag(TT,eye(t2));
      struc=[struc;st(kk,:)];
      Jeig=[Jeig;RJ(kk,:)];
   else
      T2=[];
      for kkk=1:sum(st(kk,:)~=0)
         t3=sum(st(kk,1:kkk-1));t4=st(kk,kkk);
         At=AA(t1+t3+1:t1+t3+t4,t1+t3+1:t1+t3+t4);
         [At,t5,t5,Z]=qz(At,eye(t4));
         AA(t1+t3+1:t1+t3+t4,t1+t3+1:t1+t3+t4)=At;
         T2=blkdiag(T2,Z);
      end
      T1=[];
      for kkk=1:t2/2
         T1(kkk,2*kkk-1)=1;
         T1(t2/2+kkk,2*kkk)=1;
      end
      TT=blkdiag(TT,T2*T1');
      AA(t1+1:t1+t2,t1+1:t1+t2)=T1*AA(t1+1:t1+t2,t1+1:t1+t2)*T1';
      a1=AA(t1+1,t1+1);a2=AA(t1+t2,t1+t2);
      Jeig=[Jeig;[real(a1),imag(a1)];[real(a2),imag(a2)]];
      struc=[struc;st(kk,:)/2;st(kk,:)/2];
   end
end

if size(TT,1)~=n
   AA=A;T=eye(n);struc=[];Jeig=[];err_of_JCF=inf;
   return;
end
T=T*TT;

t1=eye(n)+diag(ones(1,n-1),1);
AA=AA.*t1;
err_of_JCF=norm(inv(T)*A*T-AA);