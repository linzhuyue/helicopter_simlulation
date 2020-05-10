function gam=gm8star(A,B,C,D,E,tol)

% gms8 = gm8star(A,B,C2,D2,E[,tol])
%
% Calculates the infimum or the best achievable performance of the H_infty 
% control problem for the system:
%   .
%   x =  A x +  B u + E w
%   y =  x
%   h = C2 x + D2 u
%
% under all possible stabilizing state feedback controllers.
%
% See also gm8sos, h8state, h8out, gm2star and dgm2star.

% Modified by Ben M. Chen on April 17, 2020 at CUHK

if nargin==5
   tol=1e-8;
end

dc=0;
[A,B,C,D,Gs,Go,Gi,dims,lv,rv,qv,m0,tt,for_gm8_Gmor]=scb(A,B,C,D,tol,dc);

nan=dims(1);na0=dims(2);nap=dims(3);na=sum(dims(1:3));
nb=dims(4);nc=dims(5);nd=dims(6);
n=sum(dims);
pb=length(lv);mc=length(rv);md=length(qv);
[p,m]=size(D);

if na0~=0
   disp(' ')
   disp('...There are invariant zeros on the imaginary axis...')
   disp('...Try to use another m-function gm8sos.m instead...')
   gam=Inf;
   return;
end

if (nap+nb)==0
   gam=0;
   return;
end

Ass=A(nan+na0+1:na+nb,nan+na0+1:na+nb);
B0s=B(nan+na0+1:na+nb,1:m0);
Cd=C(m0+1:m0+md,n-nd+1:n);
Lsd=A(nan+na0+1:na+nb,n-nd+1:n)*Cd';
Bs=[B0s Lsd];

Cs=[zeros(m0+md,nap+nb);zeros(pb,nap),C(m0+md+1:p,na+1:na+nb)];
Cs=Go*Cs;
Ds=Go*[eye(m0+md);zeros(p-m0-md,m0+md)];

tE=inv(Gs)*E;
Es=tE(nan+na0+1:na+nb,:);

dt=inv(Ds'*Ds);
At=Ass-Bs*dt*Ds'*Cs;
bt=Bs*dt*Bs';
et=Es*Es';
Ct=Cs'*Cs-Cs'*Ds*dt*Ds'*Cs;

% Directly compute Infima of H-infinity optimization if possible ......
C21=for_gm8_Gmor*C(m0+1:p,nan+na0+1:na+nb);
C23=for_gm8_Gmor*[Cd*Cd';zeros(pb,md)];
C23i=inv(C23'*C23);
Ax=Ass-Lsd*C23i*C23'*C21;
BBx=B0s*B0s'+Lsd*C23i*Lsd';
CCx=C21'*C21-C21'*C23*C23i*C23'*C21;

Aaap=A(nan+na0+1:na,nan+na0+1:na);
Eap=Es(1:nap,:);
invSx=are(Ax,BBx,CCx);
if all(all(invSx~=inf))
   Tax=lyap(Aaap,-Eap*Eap');
   Tx=blkdiag(Tax,zeros(nb,nb));
   n_l=sqrt(max(eig(Tx*invSx)));
   if norm(Es(nap+1:nap+nb,:))<tol
      gam=n_l;
      return;
   end
end

n_l=0;n_h=n_l+1;
bew=0;
while (n_h-n_l)>20000*tol*n_l
   if bew==0
      n_h=10*n_h;
   end
   gam=(n_l+n_h)/2;
   Bt=bt-et/gam/gam;
   [X,flag]=zzare(At,Bt,Ct);
   if flag==1
      n_h=gam;
      bew=1;
   else
      n_l=gam;
   end
end
