function [AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0,err_scb,for_gm8_Gmor]=zzscb(A,B,C,D,tol,dc)

%ZZSCB  Special Coordinate Basis for Continuous-time Systems
%
%     [At,Bt,Ct,Dt,Gms,Gmo,Gmi,dim] = zzscb(A,B,C,D)
%
%     decomposes a continuous-time system characterized by (A,B,C,D)
%     into the standard SCB form with state subspaces x_a being
%     separated into stable, marginally stable and unstable parts
%     (in continuous-time sense), and x_d being decomposed into the
%     form of integration chains.
%
%     Input Parameters:
%           .
%           x = A x + B u,   y = C x + D u
%
%     Output Parameters:
%           .
%           x_t = At x_t + Bt u_t,   y_t = Ct x_t + Dt u_t
%
%     where x_t = [ x_a^-  x_a^0  x_a^+  x_b  x_c  x_d ]' with
%
%           dim = [ n_a^-, n_a^0, n_a^+, n_b, n_c, n_d ],
%
%     and Gms, Gmo & Gmi = state, output & input transformations.
%
%     See also SCBRAW, DSCB and SSD.

%     dc : dc=0, for continuous-time system; (default)
%          dc=1, for discrete-time system.

if nargin==4,
   tol=1e-8; dc=0;
elseif nargin==5
   dc=0;
end

n=size(A,1);
pp=size(C,1);
mm=size(B,2);

%Check whether [B;D] and [C,D] are full rank.
[u0,s0,v0]=svd([B;D]);
m=rank(s0,tol);
if mm==m,
   v0=eye(m);
else %[B;D] is not full rank
   B=B*v0;B=B(:,1:m);
   if ~isempty(D)
      D=D*v0;D=D(:,1:m);
   end
end;
[u1,s1,v1]=svd([C,D]);
p=rank(s1,tol);
if pp==p,
   u1=eye(p);
else%*** [C,D] is not full rank
   C=u1'*C;C=C(1:p,:);
   if ~isempty(D)
      D=u1'*D;D=D(1:p,:);
   end
end;

%Compute m0
if norm(D)>tol
   m0=rank(D,tol);
   [u,s,v]=svd(D); tt=[]; ttt=[];
   for kk=1:m0,
      tt(kk,kk)=sqrt(s(kk,kk));
      ttt(kk,kk)=1/sqrt(s(kk,kk));
      tt2(kk,kk)=s(kk,kk);
      ttt2(kk,kk)=1/s(kk,kk);
   end
   T2=u*blkdiag(tt,eye(p-m0));
   invT2=blkdiag(ttt,eye(p-m0))*u';
   T3=v*blkdiag(ttt,eye(m-m0));
   invT3=v'*blkdiag(tt,eye(m-m0));
   if nargin==7
      T2=u;
      invT2=u';
      T3=v*blkdiag(ttt2,eye(m-m0));
      invT3=v'*blkdiag(tt2,eye(m-m0));
   end
   Bt=B*T3;
   B0=Bt(:,1:m0);
   B1=Bt(:,m0+1:m);
   Ct=invT2*C;
   C0=Ct(1:m0,:);
   C1=Ct(m0+1:p,:);
   A1=A-B0*C0;
elseif norm(D)<=tol
   A1=A;B1=B;C1=C;
   B0=zeros(n,0);C0=zeros(0,n);D0=[];
   m0=0;T2=eye(p);invT2=eye(p);T3=eye(m);invT3=eye(m);
end

D0=eye(m0);

%Structure decomposion of (A1,B1,C1)
rv=[]; lv=[]; qv=[]; Go=[]; Gi=[];
if (isempty(C1)) & (isempty(B1))
   AA=A1; BB=[]; CC=[]; di3=[0,0,0]; Gs=eye(size(A1,1));
end
if (~isempty(C1)) & (isempty(B1))
   [AA,CC,Gs,Go,no,lv,k,k,invGs]=osd(A1,C1,tol);
   BB=[];di3=[sum(lv),0,0];
end
if (isempty(C1)) & (~isempty(B1))
   [AA,BB,Gs,Gi,no,rv,k,k,invGs]=csd(A1,B1,tol);
   di3=[0,sum(rv),0];CC=[];
end
if (isempty(B1) | isempty(C1))
   nt=size(AA,1)-sum(di3);
   At=AA(1:nt,1:nt);
   if dc==0
      [Aa,tt,nn,n0,np]=ssd(At,tol);
   else
      [Aa,tt,nn,n0,np]=dssd(At,tol);
   end
   t1=blkdiag(tt,eye(n-nt));
   Gs=Gs*t1; invGs=inv(Gs);
   AA(1:nt,1:nt)=Aa; dims=[nn,n0,np,di3];
end
if (~isempty(C1)) & (~isempty(B1))
   [AA,BB,CC,Gs,Go,Gi,lv,rv,qv,dims,invGs]=zzspscb(A1,B1,C1,tol,dc);
end
DD=zeros(size(CC,1),size(BB,2));

na=sum(dims(1:3));
AA(1:na,1:na)=AA(1:na,1:na).*blkdiag(ones(dims(1)),ones(dims(2)),ones(dims(3)));

%Structure decomposition of (A,B,C,D)
B0=invGs*B0;C0=C0*Gs;
BB=[B0,BB]; CC=[C0;CC];
DD=zeros(p,m); DD(1:m0,1:m0)=eye(m0);
for_gm8_Gmor=Go;
Go=T2*blkdiag(eye(m0),Go);
Gi=T3*blkdiag(eye(m0),Gi);

err_scb_A=norm(AA+B0*C0-invGs*A*Gs);err_scb_B=0;err_scb_C=0;
if ~isempty(B)
   err_scb_B=norm(BB-invGs*B*Gi);
end
if ~isempty(C)
   err_scb_C=norm(CC-inv(Go)*C*Gs);
end

%If [B;D] and [C,D] are not of full rank.
BB=[BB,zeros(n,mm-m)];
CC=[CC;zeros(pp-p,n)];
DD=blkdiag(DD,zeros(pp-p,mm-m));
Go=u1*blkdiag(Go,eye(pp-p));
Gi=v0*blkdiag(Gi,eye(mm-m));
err_scb=[err_scb_A,err_scb_B,err_scb_C];