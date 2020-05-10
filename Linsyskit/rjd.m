function [J,T,stru,Jeig,err_of_RJD,tol,flag]=rjd(A,tol,ordd)

%RJD  Real Jordan Decomposition
%
%     [J,T] = rjd(A[,tol])
%
%     generates a transformation that transforms a square matrix into
%     the real Jordan canonical form, i.e.,
%
%                           [ J_1         ]
%         inv(T)*A*T = J =  |     .       |
%                           |       .     |
%                           [         J_k ]
%
%      where each block J_i, i=1, ..., k, has the following form:
%
%                [ eig_i    1          ]    [ Eig_i    1          ]
%                |     .       .       |    |     .       .       |
%         J_i =  |       .       .     | or |       .       .     |
%                |         eig_i   1   |    |         Eig_i   1   |
%                [               eig_i ]    [               Eig_i ]
%
%      for real eig_i or for eig_i = mu_i + j * omiga_i, for which
%
%                Eig_i =  [   mu_i   omiga_i ]
%                         [ -omiga_i   mu_i  ]
%
%      See also JCF.

%     Note that: if ordd=1, the sizes of jordan blocks with the same eigenvalue will be in decending order.
%     stru(i,:)= the sizes of jordan blocks with the eigenvalue (mu_i + j * omega_i) and (mu_i - j * omega_i);
%     Jeig(i,:)=[mu_i, omega_i];

%     flag=1, the decomposition is reasonable for the given tol.
%     flag=0, the decomposition can not be carried out for the given tol.

if nargin==1
   tol=1234567;ordd=0;
elseif nargin==2
   ordd=0;
end

%if nargin~=1
if tol~=1234567
   [J,T,stru,Jeig,err_of_RJD,flag]=zzrjd(A,tol,ordd);
   return
else
   J=A;T=eye(size(A,1));stru=[];Jeig=zeros(0,2);err_of_RJD=inf;tol=1;flag=0;
   for k=1:8
      to=10^(-k);
      [J1,T1,stru1,Jeig1,err_of_RJD1,flag1]=zzrjd(A,to,ordd);
      if err_of_RJD>err_of_RJD1
         J=J1;
         T=T1;
         stru=stru1;
         Jeig=Jeig1;
         err_of_RJD=err_of_RJD1;
         tol=to;
         flag=flag1;
      end
   end
end