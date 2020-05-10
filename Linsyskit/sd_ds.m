function [Es4,As4,Bs4,Cs4,Ds4,Ez4,B0,C0,L4,Ck4,Dk4,yrb4,Ge,Gs,Go,invGi,dims,lv,rv,qv,na,m0]=sd_ds(E,A,B,C,D,tol)

%SD_DS  Structural Decomposition of Continuous-time Descriptor System
%
%       [Es,As0,Bs,Cs,Ds,Ez,B0,C0,Psi,Psc,Psd,Psr,Gme,Gms,Gmo,iGmi,dim]= sd_ds(E,A,B,C,D)
%
%       generates the structural decomposition of a descriptor system
%       characterized by (E,A,B,C,D).
%
%       Input Parameters:
%             .
%           E x = A x + B u,   y = C x + D u
%
%       Output Parameters:
%             .
%          Es x_t = As x_t + Bs u_t,   y_t = Cs x_t + Ds u_t
%
%       where As = As0 + B0 C0 and x_t = [ x_z  x_e  x_a  x_b  x_c  x_d ]' with
%
%             dim = [ n_z, n_e, n_a, n_b, n_c, n_d ],
%
%       (Es,As,Bs,Cs,Ds) has the same transfer function as that of
%       the original system. Ez, Psi, Psc, Psd and Psr are some
%       matrices or vectors whose elements are either polynomials
%       or rational functions of s. In particular,
%
%       Psc * x_t + Psd * u_t = Psr * x_z.
%
%       Gms, Gmo & iGmi = state, output transformations and the
%       inverse of the input transformation, and finally, Gme is
%       a nonsingular transformation on matrix E. Both Gme and
%       iGmi have their entries being some polynomials of s.
%
%       Note: Applicable to discrete-time descriptor systems too.
%
%       See also EA_SD, SCB and DSCB.

if nargin==5
   tol=1e-7;
end

%STEP MIMO-SDDS.1
n=size(A,1);
m=size(B,2);
p=size(C,1);

[PEQ,PAQ,P,Q,n1,n2]=ea_ds(E,A,tol);

N=PEQ(n1+1:n,n1+1:n);A1=PAQ(1:n1,1:n1);
PB=P*B;
CQ=C*Q;
if ~isempty(B)
   B1=PB(1:n1,:);
   B2=PB(n1+1:n,:);
else
   B1=zeros(n1,0);B2=zeros(n-n1,0);
end
if ~isempty(C)
   C1=CQ(:,1:n1);
   C2=CQ(:,n1+1:n);
else
   C1=zeros(0,n1);C2=zeros(0,n-n1);
end

%STEP MIMO-SDDS.2
[Ad,Bd,Cd,t2,t1]=ctrbf(N,B2,zeros(0,n2),tol);
t2=inv(t2);
nv=sum(t1);nz=n2-nv;
if nz==0
   Ad=N;
   Bd=B2;
   t2=eye(n2);
end

[At,Bt,tv,Ti,zpi]=bdcsd(Ad(nz+1:n2,nz+1:n2),Bd(nz+1:n2,:),tol,7654321);
if zpi==0
   zpi=[];
end
%[At,Bt,tv,Ti,zpi]=zzbdcsdmax(Ad(nz+1:n2,nz+1:n2),Bd(nz+1:n2,:),tol);
%N,B2,Ad,Bd,nz,break

ne=length(zpi);
nt=sum(zpi>1);
for k=1:ne
    wi(k)=sum(zpi(1:k));
    w2(k)=n1+sum(zpi(1:k-1));
    w3(k)=n1+sum(zpi(1:k));
end

t3=blkdiag(eye(nz),tv);
t4=eye(n2);t4=[t4(nz+1:n2,:);t4(1:nz,:)]';
Ts=t2*t3*t4;

Btt=inv(Ts)*B2*Ti;
t1=eye(m);

for k=1:ne
   t1(k,:)=Btt(wi(k),:);
   At(wi(k),:)=0;
   Bt(wi(k),:)=0;Bt(wi(k),k)=1;
end
t1=inv(t1);Ti=Ti*t1;
hN=[[At;zeros(nz,nv)],[inv(tv)*Ad(nz+1:n2,1:nz);Ad(1:nz,1:nz)]];
hB2=[Bt;zeros(nz,m)];

Tt=blkdiag(eye(n1),Ts);
Ge0=inv(Tt)*P;
Gs0=Q*Tt;
Gi0=Ti;

Et0=blkdiag(eye(n1),hN);
At0=blkdiag(A1,eye(n2));
Bt0=[B1*Ti;hB2];
Ct0=[C1,C2*Ts];
Dt0=D*Ti;

%Et0,At0,Bt0,Ct0,Dt0,break

s=sym('s');
t=sym(eye(n));
Ge1=sym(eye(n));

for k=1:nt
   if w2(k)>0
      Ge1(1:w2(k),w3(k))=-Bt0(1:w2(k),k);
   end
   for kk=1:zpi(k)
      Ge1(w3(k),w2(k)+kk)=s^(kk-1);
   end
end
Et1=Ge1*Et0;
At1=Ge1*At0;
Bt1=Ge1*Bt0;
Psiall=sym(zeros(ne,m));
Bs1=Bt1;
for k=1:ne
   Psiall(k,:)=Bt1(w3(k),:);
   Bs1(w3(k),:)=0; Bs1(w3(k),k)=1;
end

if isempty(Bs1)
   Bs1=[];
else
   Bs1=subs(Bs1,0);
end

invGi1=[Psiall;zeros(m-ne,ne),eye(m-ne)];

L1=sym(zeros(n,n));Ez1=sym(zeros(n,n));
for k=1:nt
   L1(w3(k),w2(k)+1:w3(k))=Et1(w3(k),w2(k)+1:w3(k));
end
Ez1(1:n-nz,n-nz+1:n)=Et1(1:n-nz,n-nz+1:n);
Es1=Et1-L1-Ez1;
As1=At1-s*L1;

if isempty(Es1)
   Es1=[];
else
   Es1=subs(Es1,0);
end

if isempty(As1)
   As1=[];
else
   As1=subs(As1,0);
end

Ge2=eye(n);
for k=1:nt
   Ge2(w2(k)+1,w3(k))=-1;
end

Es2=Ge2*Es1;
Ez2=Ge2*Ez1;
L2=Ge2*L1;
As2=Ge2*As1;
Bs2=Ge2*Bs1;
%Es2,Ez2,L2,As2,Bs2,break

%For Output
Cs2=Ct0;Ck2=zeros(size(Cs2));
for k=1:nt
   Cs2(:,w2(k)+1)=0;
   Ck2(:,w2(k)+1)=Ct0(:,w2(k)+1);
   Cs2(:,w3(k))=Ct0(:,w3(k))-Dt0(:,k);
   Ck2(:,w3(k))=Dt0(:,k);
end

t1=n1+sum(zpi(1:nt));
t2=n1+sum(zpi(1:ne));
Cs2(:,t1+1:t2)=0;
Ck2(:,t1+1:t2)=Ct0(:,t1+1:t2);

Ds2=Dt0;Dk2=zeros(size(Dt0));Dk2t=zeros(p,ne);
for k=1:nt
   Ds2(:,k)=-Ct0(:,w2(k)+1);
   Dk2(:,k)=Ct0(:,w2(k)+1);
end
for k=nt+1:ne
   Ds2(:,k)=Dt0(:,k)-Ct0(:,w2(k)+1);
   Dk2(:,k)=Ct0(:,w2(k)+1);
end
Dk2t=[Dt0(:,1:nt),zeros(p,ne-nt)];
if ne~=0
   if m~=ne
      Dk2=Dk2+Dk2t*inv(Psiall(:,1:ne))*[eye(ne),-Psiall(:,ne+1:m)];
   else
      Dk2=Dk2+Dk2t*inv(Psiall(:,1:ne));
   end
end

Ckr1=[Dt0(:,1:nt),zeros(p,ne-nt)];
Ckr2=zeros(p,ne);
for k=1:ne
   Ckr2(:,k)=Ct0(:,w2(k)+1);
end

yrb2=sym(zeros(p,nz));
if nz>0
   sr1=sym(zeros(ne,nz));sr2=sym(zeros(ne,nz));
   for k=1:ne
      sr1(k,:)=s*Et0(w3(k),n-nz+1:n);
      sr2(k,:)=s*Ez1(w3(k),n-nz+1:n);
   end
   yrb2=Ckr1*sr1+Ckr2*sr2;
end

t=eye(n);
Ge3=t(n-nz+1:n,:);
Gs3=t(n-nz+1:n,:);
for k=1:ne
    Ge3=[Ge3;t(n1+wi(k),:)];
    Gs3=[Gs3;t(n1+1+sum(zpi(1:k-1)),:)];
end
Ge3=[Ge3;t(1:n1,:)];
Gs3=[Gs3;t(1:n1,:)];
for k=1:nt
    Ge3=[Ge3;t(sum(zpi(1:k-1))+n1+1:n1+wi(k)-1,:)];
    Gs3=[Gs3;t(n1+2+sum(zpi(1:k-1)):n1+wi(k),:)];
end
Gs3=Gs3';

Es3=Ge3*Es2*Gs3;
As3=Ge3*As2*Gs3;
Bs3=Ge3*Bs2;
Cs3=Cs2*Gs3;
Ds3=Ds2;

Ez3=Ge3*Ez2*Gs3;
L3=Ge3*L2*Gs3;
Ck3=Ck2*Gs3;
Dk3=Dk2;
yrb3=yrb2;

%STEP MIMO-SDDS.3
bA=As3(nz+ne+1:n,nz+ne+1:n);
if isempty(Bs3)
   bB=[];
else
   bB=Bs3(nz+ne+1:n,:);
end
if isempty(Cs3)
   bC=[];
else
   bC=Cs3(:,nz+ne+1:n);
end
bD=Ds3;

dc=0;
[At,Bt,Ct,Dt,Tt,Go4,Gi4,dims,lv,rv,qv,m0]=scb(bA,bB,bC,bD,tol);

na=dims(1:3);
dims=[nz,ne,sum(na),dims(4:6)];

Ge4=blkdiag(eye(nz+ne),inv(Tt));
Gs4=blkdiag(eye(nz+ne),Tt);

Es4=Es3;
As4=Ge4*As3*Gs4;
Bs4=Ge4*Bs3*Gi4;
Cs4=inv(Go4)*Cs3*Gs4;
Ds4=inv(Go4)*Ds3*Gi4;
As4(nz+ne+1:n,nz+ne+1:n)=At;
Bs4(nz+ne+1:n,:)=Bt;
Cs4(:,nz+ne+1:n)=Ct;
Ds4=Dt;

Ez4=Ge4*Ez3*Gs4;
L4=Ge4*L3*Gs4;
Ck4=inv(Go4)*Ck3*Gs4;
Dk4=inv(Go4)*Dk3*Gi4;
yrb4=[];
if size(yrb3,1)~=0
   yrb4=inv(Go4)*yrb3;
end

Ge=Ge4*Ge3*Ge2*Ge1*Ge0;
Gs=Gs0*Gs3*Gs4;
Gi=Ti*inv(invGi1)*Gi4;
invGi=inv(Gi4)*invGi1*inv(Ti);
Go=Go4;

dig=10;
Ez4=vpa(Ez4,dig);
L4=vpa(L4,dig);

Dk4=vpa(Dk4,dig);
yrb4=vpa(yrb4,dig);
Ge=vpa(Ge,dig);
invGi=vpa(invGi,dig);

B0=Bs4(:,1:m0);B0(1:ne+nz,:)=0;
C0=Cs4(1:m0,:);C0(:,1:ne+nz)=0;

%subs(Ge*E*Gs-Es4-Ez4-L4,2)
%subs(Ge*A*Gs-As4-B0*C0-s*L4,2)
%subs(Ge*B*Gi-Bs4,2)
%subs(inv(Go)*C*Gs-Cs4-Ck4,2)
%subs(inv(Go)*D*Gi-Ds4-Dk4,2)
