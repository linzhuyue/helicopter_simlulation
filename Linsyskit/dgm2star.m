function gms2 = dgm2star(A,B,C,D,E,tol)

%%DGM2STAR  Infimum or Optimal Value for Discrete-time H2 Control
%
%          gms2 = dgm2star(A,B,C,D,E)
%
%          calculates the infimum or the best achievable performance
%          of the H2 suboptimal control problem for the system:
%
%              x(k+1) = A x(k) + B u(k) + E w(k)
%               h(k)  = C x(k) + D u(k)
%
%          under all possible stabilizing state feedback controllers.
%
%          See also DH2STATE, DGM8STAR and GM2STAR.

if nargin==5
   tol=1e-8;
end
dc=1;
[A,B,C,D,Gs,Go,Gi,dims,lv,rv,qv,m0]=scb(A,B,C,D,tol,dc);

%type=2;dc=1;d11_eye=1;
%[A,B,C,D,Gs,Go,Gi,dims,lv,rv,qv,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);

nn=dims(1);no=dims(2);np=dims(3);na=nn+no+np;
nb=dims(4);nc=dims(5);nd=dims(6);n=sum(dims);
pb=length(lv);mc=length(rv);md=length(qv);
[p,m]=size(D);
Bd=B(n-nd+1:n,m0+1:m0+md);
Bc=B(na+nb+1:na+nb+nc,m0+md+1:m);

if no~=0,
   disp('Sorry. There are invariant zeros on the unit circle.')
   gms2=inf;
   return;
end

Gt=eye(n);
Gt=[Gt(1:nn,:);Gt(na+nb+1:na+nb+nc,:);Gt(nn+1:na+nb,:);Gt(n-nd+1:n,:)];
Gs=Gs*Gt';
A=Gt*A*Gt';
B=Gt*B;
C=C*Gt';

Ass=A(nn+nc+1:n,nn+nc+1:n);
Bs=B(nn+nc+1:n,1:m0+md);
Cs=C(:,nn+nc+1:n);Cs(1:m0,:)=0;
Cs=Go*Cs;
Ds=zeros(p,m0+md);Ds(1:m0,1:m0)=eye(m0);
Ds=Go*Ds;

tE=inv(Gs)*E;
Es=tE(nn+nc+1:n,:);

Ps=h2dare(Ass,Bs,Cs,Ds);
%Fs=inv(Ds'*Ds+Bs'*Ps*Bs)*(Bs'*Ps*Ass+Ds'*Cs);

gms2=sqrt(trace(Es'*Ps*Es));
