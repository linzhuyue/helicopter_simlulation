function [EE,AA,U,V,nn,ll,na0]=zzfiguptri(EE,AA,tol,type);

%Generalized upper triangular form of the pair (E, A)
%
%     [EE,AA,U,V,nn,ll,na0]=zzfiguptri(E,A,tol);
%
%                    ll(1)     ll(2)       na0    ll(3)-na0     ll(4)     ll(5)
%                  [sE11-A11  sE12-A12  sE13-A13  sE14-A14  sE15-A15  sE16-A16]   nn(1)
%                  [    0     sE22-A22  sE23-A23  sE24-A24  sE25-A25  sE26-A26]   nn(2)
%     U*(sE-A)*V = [    0         0     sE33-A33  sE34-A34  sE35-A35  sE36-A36]   na0
%                  [    0         0         0     sE44-A44  sE45-A45  sE46-A46]   nn(3)-na0
%                  [    0         0         0         0     sE55-A55      0   ]   nn(4)
%                  [    0         0         0         0     sE65-A65  sE66-A66]   nn(5)
%
%     where E11 is of full-row rank;
%           E22, E33 and E44 are nonsigular;
%           E55 is of full-column rank;
%           sE11-A11 is of full-row rank;
%           eig(inv(E22)*A22)<=-tol;
%           -tol<eig(inv(E33)*A33)<=tol;
%           eig(inv(E44)*A44)=>tol;
%           sE55-A55 is of full-column rank;
%           sE66-A66 is of nonsingular;
%

if nargin==2
   tol=1e-8;type=0;
elseif nargin==3
   type=0;
end   

%E=EE;A=AA;

[EE,AA,U1,V1,nn,ll]=zzguptri2(EE,AA,tol);
Et=EE(nn(1)+1:nn(1)+nn(2),ll(1)+1:ll(1)+ll(2));
At=AA(nn(1)+1:nn(1)+nn(2),ll(1)+1:ll(1)+ll(2));

%Decomposition of Aa
na0=0;k=size(Et,1);
[Et,At,U2,V2,mu]=zzschur2new(Et,At,tol);
if (type==1)&((k-mu)>0)
   Ett=Et(mu+1:k,mu+1:k);
   Att=-At(mu+1:k,mu+1:k);
   [Ett,Att,U2t,V2t,mup]=zzschur2new(Ett,Att,tol);
   na0=k-mu-mup;
   ut=eye(k-mu);ut=[ut(mup+1:na0+mup,:);ut(1:mup,:)];
   U2t=ut*U2t;V2t=V2t*ut';
   U2=blkdiag(eye(mu),U2t)*U2;
   V2=V2*blkdiag(eye(mu),V2t);
end

Ut=eye(nn(1)+nn(2));Vt=eye(ll(1)+ll(2));
Ut(nn(1)+nn(2)+1:nn(1)+nn(2)+nn(4),nn(1)+nn(2)+nn(3)+1:sum(nn))=eye(nn(4));
Ut(nn(1)+nn(2)+nn(4)+1:sum(nn),nn(1)+nn(2)+1:nn(1)+nn(2)+nn(3))=eye(nn(3));
Vt(ll(1)+ll(2)+1:ll(1)+ll(2)+ll(3),ll(1)+ll(2)+ll(4)+1:sum(ll))=eye(ll(3));
Vt(ll(1)+ll(2)+ll(3)+1:sum(ll),ll(1)+ll(2)+1:ll(1)+ll(2)+ll(4))=eye(ll(4));
Ut(nn(1)+1:nn(1)+nn(2),nn(1)+1:nn(1)+nn(2))=U2;
Vt(ll(1)+1:ll(1)+ll(2),ll(1)+1:ll(1)+ll(2))=V2;
EE=Ut*EE*Vt;AA=Ut*AA*Vt;
EE(nn(1)+1:nn(1)+nn(2),ll(1)+1:ll(1)+ll(2))=Et;
AA(nn(1)+1:nn(1)+nn(2),ll(1)+1:ll(1)+ll(2))=At;
U=Ut*U1;V=V1*Vt;
nn=[nn(1),mu,nn(2)-mu,nn(4),nn(3)];
ll=[ll(1),mu,ll(2)-mu,ll(4),ll(3)];


%[n,ell]=size(A);
%At=ones(n,n);
%At(nn(1)+1:n,1:ll(1))=0;
%At(nn(1)+nn(2)+1:n,ll(1)+1:ll(1)+ll(2))=0;
%At(nn(1)+nn(2)+nn(3)+1:n,ll(1)+ll(2)+1:ll(1)+ll(2)+ll(3))=0;
%At(nn(1)+nn(2)+nn(3)+1:n-nn(5),ell-ll(5):ell)=0;

%err_figuptri_E=norm(U*E*V.*At-EE);
%err_figuptri_A=norm(U*A*V.*At-AA);
%err_firguptri=[err_figuptri_E,err_figuptri_A]