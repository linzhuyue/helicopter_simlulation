function [AA,T,struc,Jeig,err_of_JCF,tol,flag]=jcf(A,tol)

%JCF  Jordan Canonical Form
%
%     [J,T] = jcf(A[,tol])
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

if nargin==2
   [AA,T,struc,Jeig,err_of_JCF,flag]=zzjcf(A,tol);
   return
else
   AA=A;T=eye(size(A,1));struc=[];Jeig=zeros(0,2);err_of_JCF=inf;tol=1;flag=0;
   for k=1:8
      to=10^(-k);
      [AA1,T1,struc1,Jeig1,err_of_JCF1,flag1]=zzjcf(A,to);
      if err_of_JCF>err_of_JCF1
         AA=AA1;
         T=T1;
         struc=struc1;
         Jeig=Jeig1;
         err_of_JCF=err_of_JCF1;
         tol=to;
         flag=flag1;
      end
   end
end
