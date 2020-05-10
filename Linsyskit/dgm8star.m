function gam = dgm8star(A,B,C,D,E)

%DGM8STAR  Infimum or Optimal Value for Discrete H-infinity Control
%
%          gms8 = dgm8star(A,B,C,D,E)
%
%          calculates the infimum or the best achievable performance
%          of H-infinity suboptimal control problem for the plant:
%
%              x(k+1) = A x(k) + B u(k) + E w(k)
%               h(k)  = C x(k) + D u(k)
%
%          under all possible stabilizing state feedback controllers.
%
%          See also DH8ARE, DH8STATE, DGM2STAR and GM8STAR.

if nargin==5
   tol=1e-8;
end
dc=1;
[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=scb(A,B,C,D,tol,dc);

nn=dims(1);no=dims(2);np=dims(3);na=nn+no+np;
nb=dims(4);nc=dims(5);nd=dims(6);n=sum(dims);
pb=length(lv);mc=length(rv);md=length(qv);
p=size(C,1);
m=size(B,2);
Bd=BB(n-nd+1:n,m0+1:m0+md);
Bc=BB(na+nb+1:na+nb+nc,m0+md+1:m);

gam=[];
if no~=0,
   disp('......There are invariant zeros on the unit circle......')
   return;
end

Gt=eye(n);
Gt=[Gt(1:nn,:);Gt(na+nb+1:na+nb+nc,:);Gt(nn+1:na+nb,:);Gt(n-nd+1:n,:)];
Gs=Gs*Gt';
AA=Gt*AA*Gt';
BB=Gt*BB;
CC=CC*Gt';

Ass=AA(nn+nc+1:n,nn+nc+1:n);
Bs=BB(nn+nc+1:n,1:m0+md);
Cs=CC(:,nn+nc+1:n);Cs(1:m0,:)=0;
Cs=Go*Cs;
Ds=zeros(p,m0+md);Ds(1:m0,1:m0)=eye(m0);
Ds=Go*Ds;

tE=inv(Gs)*E;
Es=tE(nn+nc+1:n,:);

Q=E*E';

n_l=0; n_h=12345;
bew=0;

if max(size(Ass))==0
   gam=0;
else while (n_h-n_l)>2*tol*n_l
       if bew==0
          n_h=10*n_h;
       end
       gam=(n_l+n_h)/2;
       [Ps,flag]=h8dare(Ass,Bs,Cs,Ds,Es,gam);
       if flag==1
          P=inv(Gs)'*[zeros(nn+nc,n);zeros(n-nn-nc,nn+nc) Ps]*inv(Gs);      
          if max(eig(P*Q))<gam*gam
             n_h=gam;
             bew=1;
          else
             n_l=gam;
          end
       else
          n_l=gam;
       end
     end
 end