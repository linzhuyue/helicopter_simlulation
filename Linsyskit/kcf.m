function [UPV,U,V,dims,UPV_1,UPV_2]=kcf(A,B,C,D,tol)

%KCF  Kronecker Canonical Form of System Matrices
% 
%     [Ks,U,V,dims] = kcf(A,B,C,D)
% 
%     returns the Kronecker canonical form of a system matrix for
%     a system characterized by (A,B,C,D), i.e.,
% 
%                 [ sI-A -B ]
%         Ks =  U |         | V
%                 [ C     D ]
% 
%     where Ks is in the Kronecker canonical form.
% 
%                     [ na        na ]
%                     | nb+pb     nb |
%              dims = | nc     nc+mc |
%                     | nd+md  nd+md |
%                     | m0        m0 |
%                     [ r0        c0 ]
% 
%     contains the dimensions of Kronecker blocks associated with
%     x_a, x_b, x_c, x_d, I_{m_0} and 0.
% 
%     See also SCB, DSCB and MORSEIDX.


if nargin==4
   tol=1e-8;
end

[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=scb(A,B,C,D,tol);

%type=2;dc=0;d11_eye=1;
%[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);

n=size(A,1);m=size(B,2);p=size(C,1);
nb=sum(lv);pb=length(lv);
nc=sum(rv);mc=length(rv);
nd=sum(qv);md=length(qv);
na=n-nb-nc-nd;
r0=p-m0-pb-md;
c0=m-m0-mc-md;
dims=[na,na; nb+pb,nb; nc,nc+mc; nd+md,nd+md; m0,m0; r0,c0];

[Aaa,t1]=jcf(AA(1:na,1:na),tol);
t1=blkdiag(t1,eye(n-na));
Gs=Gs*t1;
AA=inv(t1)*AA*t1;
BB=inv(t1)*BB;
CC=CC*t1;

Cb=CC(m0+md+1:p,na+1:na+nb);
Bc=BB(na+nb+1:na+nb+nc,m0+md+1:m);
Bd=BB(n-nd+1:n,m0+1:m0+md);
Cd=CC(m0+1:m0+md,n-nd+1:n);

C0=CC(1:m0,:);
Ed=Bd'*AA(n-nd+1:n,:);
Ec=Bc'*AA(na+nb+1:n-nd,:); Ec(:,[na+1:na+nb,n-nd+1:n])=0;
tF=-[C0;Ed;Ec];

B0=BB(:,1:m0);
Ld=(AA(:,n-nd+1:n)-[zeros(n-nd,nd);Bd*Ed(:,n-nd+1:n)])*Cd';
Lb=AA(:,na+1:na+nb)*Cb'; Lb(n-nd+1:n,:)=0;
tK=-[B0,Ld,Lb];

U1=[inv(Gs),-tK*inv(Go);zeros(p,n),inv(Go)];
V1=[Gs,zeros(n,m);Gi*tF,Gi];

Ut=eye(n+p);
U2=Ut(1:na,:);
for kk=1:pb
   t1=na+sum(lv(1:kk-1));
   t2=na+sum(lv(1:kk));
   U2(t1+kk,:)=-Ut(n+m0+md+kk,:);
   U2(t1+kk+1:t2+kk,:)=Ut(t1+1:t2,:);
end
U2(na+nb+pb+1:na+nb+pb+nc,:)=Ut(na+nb+1:na+nb+nc,:);
for kk=1:md
   t1=n-nd+sum(qv(1:kk-1));
   t2=n-nd+sum(qv(1:kk));
   U2(t1+pb+kk,:)=-Ut(n+m0+kk,:);
   U2(t1+pb+kk+1:t2+pb+kk,:)=Ut(t1+1:t2,:);
end
U2(n+pb+md+1:n+pb+md+m0,:)=Ut(n+1:n+m0,:);
U2(n+pb+md+m0+1:n+p,:)=Ut(n+pb+md+m0+1:n+p,:);

Vt=eye(n+m);
V2=Vt(:,1:na+nb);
for kk=1:mc
   t1=na+nb+sum(rv(1:kk-1));
   t2=na+nb+sum(rv(1:kk));
   V2(:,t1+kk:t2+kk-1)=Vt(:,t1+1:t2);
   V2(:,t2+kk)=Vt(:,n+m0+md+kk);
end
for kk=1:md
   t1=n-nd+sum(qv(1:kk-1));
   t2=n-nd+sum(qv(1:kk));
   V2(:,t1+mc+kk:t2+mc+kk-1)=Vt(:,t1+1:t2);
   V2(:,t2+mc+kk)=Vt(:,n+m0+kk);
end
V2(:,n+mc+md+1:n+mc+md+m0)=Vt(:,n+1:n+m0);
V2(:,n+mc+md+m0+1:n+m)=Vt(:,n+mc+md+m0+1:n+m);

tt=[];
for kk=1:md
   tt=blkdiag(tt,rot90(eye(qv(kk)+1)));
end
Ut(n-nd+pb+1:n+pb+md,n-nd+pb+1:n+pb+md)=-tt;
Vt(n-nd+mc+1:n+mc+md,n-nd+mc+1:n+mc+md)=tt;

U=Ut*U2*U1;
V=V1*V2*Vt;

UPV_1=eye(na); UPV_2=-Aaa;
for kk=1:pb
   UPV_1=blkdiag(UPV_1,[zeros(1,lv(kk));eye(lv(kk))]);
   UPV_2=blkdiag(UPV_2,-[eye(lv(kk));zeros(1,lv(kk))]);
end
for kk=1:mc
   UPV_1=blkdiag(UPV_1,[eye(rv(kk)),zeros(rv(kk),1)]);
   UPV_2=blkdiag(UPV_2,-[zeros(rv(kk),1),eye(rv(kk))]);
end
for kk=1:md
   UPV_1=blkdiag(UPV_1,-diag(ones(qv(kk),1),1));
   UPV_2=blkdiag(UPV_2,eye(qv(kk)+1));
end
UPV_1=blkdiag(UPV_1,zeros(m0,m0),zeros(p-pb-md-m0,m-mc-md-m0));
UPV_2=blkdiag(UPV_2,eye(m0,m0),zeros(p-pb-md-m0,m-mc-md-m0));

s=sym('s');
UPV=s*UPV_1+UPV_2;
dig=10;
UPV=vpa(UPV,dig);

%P=[-A,-B;C,D];
%pt1=zeros(size(P));pt1(1:n,1:n)=eye(n);pt1=Ut*U2*pt1*V2*Vt;
%pt2=U*P*V;
%err1=pt1-UPV_1, err2=pt2-UPV_2
