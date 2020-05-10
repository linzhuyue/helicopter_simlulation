function [J,T,stru,Jeig,err_of_RJD,flag]=zzrjd(A,tol,ordd)

%zzRJD  Real Jordan Decomposition
%
%     [J,T] = zzrjd(A[,tol])
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
   tol=1e-6;ordd=0;
elseif nargin==2
   ordd=0;
end

flag=1;
if isempty(A)
   J=[];T=[];stru=[];Jeig=zeros(0,2);err_of_RJD=0;
   return;
end

n=size(A,1);
[At,P,dims,err]=zzbdeig(A,tol);
if err==inf
   J=A;T=eye(n);stru=[];Jeig=zeros(0,2);err_of_RJD=inf;
   return
end

%S=eye(n);
S=[];
sturc=[];
Jeig=[];
ell=length(dims);

for k=1:ell 
   ni=dims(k);
   nt=sum(dims(1:k-1)); 
   Ai=At(nt+1:nt+ni,nt+1:nt+ni);
   lambda_i=eig(Ai,'nobalance');
   mu_i=mean(real(lambda_i));
   omega_i=mean(abs(imag(lambda_i)));
   if omega_i<tol
      [Ji,Pi,stru_i,flag]=zzrepvaluej(Ai,tol);
      if flag==0
         J=A;T=eye(n);stru=[];;Jeig=zeros(0,2);err_of_RJD=inf;
         return;
      end
      %S(nt+1:nt+ni,nt+1:nt+ni)=Pi;
      S=blkdiag(S,Pi);
      stru(k,1:length(stru_i))=stru_i;
      Jeig(k,:)=[mu_i, 0];
   else
      Zi=[Ai-mu_i*eye(ni), omega_i*eye(ni); -omega_i*eye(ni), Ai-mu_i*eye(ni)];
      Zit=inv(Zi-omega_i*eye(2*ni))*(Zi+omega_i*eye(2*ni));
      [t1,Si0,t,t,t]=ssd(Zit,tol);
      t1=inv(Si0)*Zi*Si0;
      Zi0=t1(1:ni,1:ni);
      [Ji,Si1,t1,flag]=zzrepvaluej(Zi0,tol);
      if flag==0
         J=A;T=eye(n);stru=[];Jeig=zeros(0,2);,err_of_RJD=inf;
         return;
      end
      stru_i=[];
      for k2=1:2:length(t1)
         stru_i=[stru_i,2*t1(k2)];
      end
      Si=Si0*blkdiag(Si1,eye(ni));
      tSi=[];
      for k1=1:length(stru_i)
         x=sum(stru_i(1:k1-1));
         for k2=1:stru_i(k1)/2
            tSi=[tSi, Si(1:ni,x+k2), Si(ni+1:2*ni,x+k2)];
         end
      end
      %S(nt+1:nt+ni,nt+1:nt+ni)=tSi;
      S=blkdiag(S,tSi);
      stru(k,1:length(stru_i))=stru_i;
      Jeig(k,:)=[mu_i, omega_i];
   end
end

if any(size(S)~=n)
   J=A;T=eye(n);stru=[];Jeig=zeros(0,2);err_of_RJD=inf;flag=0;
   return;
end
   
T=real(P*S);

if ordd==1
   jr=size(stru,1);Han=[];
   for k1=1:jr
      y1=eye(sum(stru(k1,:)));kk=sum(stru(k1,:)~=0);
      t3=[];
      for k2=1:kk
         t1=sum(stru(k1,1:k2-1));t2=sum(stru(k1,1:k2));
         t3=[y1(t1+1:t2,:);t3];
         ss(k1,k2)=stru(k1,kk+1-k2);
      end
      Han=blkdiag(Han,t3);
   end
   stru=ss;
   T=T*Han';
end

%Makes the elments should be zero be zero.
tt=eye(n);
[ell,sigmai]=size(stru); 
for k1=1:ell 
   t0=sum(sum(stru(1:k1-1,:))); 
   for k2=1:sum(stru(k1,:)~=0) 
      t1=t0+sum(stru(k1,1:k2-1));t2=t0+sum(stru(k1,1:k2)); 
      if Jeig(k1,2)==0 
         tt(t1+1:t2-1,t1+2:t2)=tt(t1+1:t2-1,t1+2:t2)+eye(stru(k1,k2)-1); 
      else 
         tt(t1+1:t2-2,t1+3:t2)=tt(t1+1:t2-2,t1+3:t2)+eye(stru(k1,k2)-2); 
         for k3=1:stru(k1,k2)/2
             tt(t1+2*k3-1:t1+2*k3,t1+2*k3-1:t1+2*k3)=ones(2,2); 
         end 
      end 
   end 
end

J=inv(T)*A*T.*tt;
err_of_RJD=norm(inv(T)*A*T-J);