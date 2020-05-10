function [AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0,err_scbchu,for_gm8_Gmor]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye)

%SCBCHU  Raw Decomposition of Proper System.
%
%     [AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=zzscbchu(A,B,C,D[,tol,type,dc,d11_eye])
%
%     decomposes the system (A,B,C,D) into its scb structure. The algorithm consist
%     of three steps: reduction by orthogonal transformations, reduction by generalized
%     Sylvester equations, and extraction of infinite zero structure.
%
%     Inputs:
%
%        A, B, C, D : state space matrices of a given system.
%              type : type=0, raw decomposition without intergration chains; (default)
%                     type=1, giving infinite zero structure;
%                     type=2, giving stable, unstable invariant zeros and I2,I3,I4;
%                dc : dc=0, for continuous-time system; (default)
%                     dc=1, for discrete-time system.
%           dll_eye : dll_eye=0,  D11 is not identity matrix;
%                     dll_eye=1,  D11 is identity matrix; (default)
%
%     Outputs:
%        AA, BB, CC, DD : state space matrices in an s.c.b.
%        Gs, Go, Gi : state, output and input transformations.
%        lv : = [l1, l2, ..., l_{pb}], observability index corresponding to x_b;
%        rv : = [r1, r2, ..., r_{mc}], ontrollability index corresponding to x_c;
%        qv : = [q1, q2, ..., q_{md}], infinite zeros, corresponding to x_d;
%        m0 : =  rank of D.
%        dims : = [nan, na0, nap, nb, nc, nd].
%
%     For details, see the manual and the references therein.
%
%     See also scb.

%     Refer to the following for detail
%     Chu DL et al, IEEE Trans. on Automatic Control, Vol.47, No.11

if nargin==4
   tol=1e-7; type=0; dc=0; d11_eye=1;
elseif nargin==5
   type=0; dc=0; d11_eye=1;
elseif nargin==6
   dc=0; d11_eye=1;
elseif nargin==7
   d11_eye=1;
end;

n=size(A,1);
pp=size(C,1);
mm=size(B,2);

%Check whether [B;D] and [C,D] are full rank.
[u0,s0,v0]=svd([B;D]);
m=rank(s0,tol);
%m=sum(diag(s0)>tol);
if mm==m,
   v0=eye(m);
else%[B;D] is not full rank
   B=B*v0;B=B(:,1:m);
   if ~isempty(D)
      D=D*v0;D=D(:,1:m);
   end
end;
[u1,s1,v1]=svd([C,D]);
p=rank(s1,tol);
%p=sum(diag(s1)>tol);
if pp==p,
   u1=eye(p);
else%[C,D] is not full rank
   C=u1'*C;C=C(1:p,:);
   if ~isempty(D)
      D=u1'*D;D=D(1:p,:);
   end
end;

if norm(D)>tol
   if d11_eye==0
      [Ud,Dt,m0]=zzrowup(D,tol);
      [Vd,D11]=zzrowup(Dt',eps,m0);Vd=Vd';D11=D11';
      Bt=B*Vd;
      B0=Bt(:,1:m0);
      BB1=Bt(:,m0+1:m);
      Ct=Ud*C;
      C0=Ct(1:m0,:);
      CCt=Ct(m0+1:p,:);
      [Q,D0]=zzrowdown([C0,D11]',eps,m0);Q=Q';D0=D0';
      EE1=Q(1:n,1:n);
      At=[A,B0;CCt,CCt*B0*0]*Q;
      AA1=At(1:n,1:n);
      CC1=At(n+1:n+p-m0,1:n);
      invUd=inv(Ud);
      invVd=inv(Vd);
   elseif d11_eye==1
      m0=rank(D,tol);
      [u,s,v]=svd(D); tt=[]; ttt=[];
      %m0=sum(diag(s)>tol);
      for kk=1:m0
         tt(kk,kk)=sqrt(s(kk,kk));
         ttt(kk,kk)=1/sqrt(s(kk,kk));
         tt2(kk,kk)=s(kk,kk);
         ttt2(kk,kk)=1/s(kk,kk);
      end;
      Ud=blkdiag(ttt,eye(p-m0))*u';
      invUd=u*blkdiag(tt,eye(p-m0));
      Vd=v*blkdiag(ttt,eye(m-m0));
      invVd=blkdiag(tt,eye(m-m0))*v';
      if nargin==9
         Ud=u';
         invUd=u;
         Vd=v*blkdiag(ttt2,eye(m-m0));
         invVd=blkdiag(tt2,eye(m-m0))*v';
      end
      Bt=B*Vd;
      B0=Bt(:,1:m0);
      BB1=Bt(:,m0+1:m);
      Ct=Ud*C;
      C0=Ct(1:m0,:);
      CC1=Ct(m0+1:p,:);
      EE1=eye(n);
      AA1=A-B0*C0;
      D11=eye(m0);
   end
elseif norm(D)<=tol
   EE1=eye(n);AA1=A;BB1=B;CC1=C;
   B0=zeros(n,0);C0=zeros(0,n);D11=[];
   m0=0;Ud=eye(p);Vd=eye(m);
   invUd=eye(p);invVd=eye(m);
end
[AA,BB,CC,V,Go,Gi,qv,dims,U]=zzescbedz(EE1,AA1,BB1,CC1,type,tol);

At=AA;
Gs=EE1*V;
invGs=U;
na=sum(dims(1:3));nb=dims(4);nc=dims(5);nd=dims(6);md=length(qv);
lv=NaN;rv=NaN;
if type==2
   Aa=AA(1:na,1:na);
   Ab=AA(na+1:na+nb,na+1:na+nb);
   Cb=CC(md+1:p-m0,na+1:na+nb);
   Ac=AA(na+nb+1:na+nb+nc,na+nb+1:na+nb+nc);
   Bc=BB(na+nb+1:na+nb+nc,md+1:m-m0);
   Ad=AA(n-nd+1:n,n-nd+1:n);
   Bd=BB(n-nd+1:n,1:md);
   Cd=CC(1:md,n-nd+1:n);
   Ta=eye(na);
   if dc==1
      [Aa,Ta,nn,no,np]=dssd(Aa,tol);
      dims(1:3)=[nn,no,np];
   elseif dc==0
      [Aa,Ta,nn,no,np]=ssd(Aa,tol);
      dims(1:3)=[nn,no,np];
   end
   %Ab,Cb,
   [Ab,Cb,Tbs,Tbo,no,lv,t1,t1,invTbs]=osd(Ab,Cb,tol);
   [Ac,Bc,Tcs,Tci,no,rv,t1,t1,invTcs]=csd(Ac,Bc,tol);
   TT=blkdiag(Ta,Tbs,Tcs,eye(nd));
   invTT=blkdiag(inv(Ta),invTbs,invTcs,eye(nd));
   Gs=Gs*TT;
   invGs=invTT*invGs;
   Go=Go*blkdiag(eye(md),Tbo);
   Gi=Gi*blkdiag(eye(md),Tci);
   AA=invTT*At*TT;

   AA(1:na,1:na)=Aa;
   AA(na+1:na+nb,na+1:na+nb)=Ab;
   CC(md+1:p-m0,na+1:na+nb)=Cb;
   AA(na+nb+1:na+nb+nc,na+nb+1:na+nb+nc)=Ac;
   BB(na+nb+1:na+nb+nc,md+1:m-m0)=Bc;
   AA(n-nd+1:n,n-nd+1:n)=Ad;
   BB(n-nd+1:n,1:md)=Bd;
   CC(1:md,n-nd+1:n)=Cd;
   AA(1:na,na+1:na+nb)=AA(1:na,na+1:na+nb)*Cb'*Cb;
   AA(na+nb+1:na+nb+nc,na+1:na+nb)=AA(na+nb+1:na+nb+nc,na+1:na+nb)*Cb'*Cb;
   AA(na+nb+1:na+nb+nc,1:na)=Bc*Bc'*AA(na+nb+1:na+nb+nc,1:na);
   AA(na+nb+1:na+nb+nc,1:na)=Bc*Bc'*AA(na+nb+1:na+nb+nc,1:na);
   
   AA(1:na,1:na)=AA(1:na,1:na).*blkdiag(ones(dims(1)),ones(dims(2)),ones(dims(3)));
end

B0=invGs*B0;C0=C0*Gs;
BB=[B0,BB]; CC=[C0;CC];
DD=blkdiag(D11,zeros(p-m0,m-m0));
for_gm8_Gmor=Go;
Go=invUd*blkdiag(eye(m0),Go);
Gi=Vd*blkdiag(eye(m0),Gi);

err_scbchu_A=norm(AA+B0*inv(D11)*C0-invGs*A*Gs);
err_scbchu_B=0;err_scbchu_C=0;
if ~isempty(B)
   err_scbchu_B=norm(BB-invGs*B*Gi);
end
if ~isempty(C)
   err_scbchu_C=norm(CC-inv(Go)*C*Gs);
end

%If [B;D] and [C,D] are not of full rank.
BB=[BB,zeros(n,mm-m)];
CC=[CC;zeros(pp-p,n)];
DD=blkdiag(DD,zeros(pp-p,mm-m));
Go=u1*blkdiag(Go,eye(pp-p));
Gi=v0*blkdiag(Gi,eye(mm-m));
err_scbchu=[err_scbchu_A,err_scbchu_B,err_scbchu_C];