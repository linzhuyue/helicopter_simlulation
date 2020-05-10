function [Am,Bm,Cm,Dm,Av,Bv,Cv,Dv]=gcfact(A,B,C,D,tol)

%GCFACT  Generalized Cascade Factorization of Continuous-time Systems
%
%        [Am,Bm,Cm,Dm,Au,Bu,Cu,Du] = gcfact(A,B,C,D)
%
%        generates a generalized cascade factorization for a system
%        (A,B,C,D) with transfer function matrix G(s), where [B' D']
%        and [C D] are assumed to be of full rank, and all the
%        'awkward' invariant zeros of (A,B,C,D) are detectable.
%
%        The generalized cascade factorization is given as
%
%                G(s) = Gm(s) U(s)
%
%        where
%
%                Gm(s) = Cm (sI - Am)^{-1} Bm + Dm
%
%        is of minimum-phase and left invertible, and
%
%                U(s) = Cu (sI - Au)^{-1} Bu + Du
%
%        is a stable right invertible and asymptotic all-pass, i.e.,
%
%                U(s) U'(-s) -> I   as |s| -> infinity
%
%        Note that users will be prompted to enter desired zero
%        locations to replace those 'awkward' invariant zeros.
%
%        See also MPFACT and IOFACT.

if nargin==4,
   tol=1e-8;
end;

dc=0;
[As,Bs,Cs,Ds,Gs,Go,Gi,dims,lv,rv,qv,m0]=scb(A,B,C,D,tol,dc);

nan=dims(1);na0=dims(2);nap=dims(3);
na=sum(dims(1:3));nb=dims(4);nc=dims(5);nd=dims(6);
p=size(C,1);
m=size(B,2);
Am=A;Cm=C;

n=size(A,1);
md=length(qv);

if (na+nc)==0,
   Bm=B;Dm=D;Av=[];Bv=zeros(0,m);Cv=zeros(m,0);Dv=eye(m);
   return;
end
Ax=[As(na+nb+1:na+nb+nc,na+nb+1:na+nb+nc),As(na+nb+1:na+nb+nc,1:na)];
Ax=[Ax;[zeros(na,nc),As(1:na,1:na)]];

Bx=zeros(na+nc,m);
if nc~=0,
   Bx(1:nc,(m0+md)+1:m)=Bs(na+nb+1:na+nb+nc,(m0+md)+1:m);
end
Bx=Bx*inv(Gi);
Bd=Bs(na+nb+nc+1:n,m0+1:(m0+md));
Cx=[Cs(1:m0,na+nb+1:na+nb+nc),Cs(1:m0,nan+1:na)];
Cx=[Cx;Bd'*[As(na+nb+nc+1:n,na+nb+1:na+nb+nc),As(na+nb+nc+1:n,1:na)]];
Dx=zeros((m0+md),m);
if (m0+md)~=0,
   Dx(:,1:(m0+md))=eye(m0+md);
end
Dx=Dx*inv(Gi);
G=inv(Gi)*inv(Gi');
invGm=G(1:(m0+md),1:(m0+md));
if norm(invGm-eye(m0+md))>tol*1e-4,
   invGm=sqrtm(invGm);
end

[a,b,c,T,k]=obsvf(Ax,Bx,Cx);
k=sum(k);
P1=eig(a(1:na+nc-k,1:na+nc-k));
if any(real(P1)>=0),
   disp(' ')
   disp('Warning : The given system is not detectable...');
   Am=[];Bm=[];Cm=[];Dm=[];Av=[];Bv=[];Cv=[];Dv=[];
   return
end
P2=eig(a(na+nc-k+1:na+nc,na+nc-k+1:na+nc));
Kx=zeros(na+nc,(m0+md));

if length(P2)~=0
   clc
   disp(' ');
   disp('* You have freedom to relocate the following invariant zeros to anywhere...');
   disp(' ');
   disp(P2');
   disp(' ');
   mmm=max(size(P2));
   
   ff=0;
   while ff==0
      disp(' ')
      disp('  The desired invariant zeros should be a self-conjugated set in a row')
      disp(' ')
      disp(['   vector with ', int2str(mmm), ' entries...'])
      disp(' ')   
      P3=input('  Enter the desired location(s) = ');
      ff=1;
      if size(P3,1)~=1 | size(P3,2)~=mmm | abs(imag(sum(P3)))>tol
         disp(' ')
         disp('  > > > The entered set is invalid. Please re-enter a valid one. < < <')
         ff=0;
      end
   end
   Kx=place(a(na+nc-k+1:na+nc,na+nc-k+1:na+nc)',c(:,na+nc-k+1:na+nc)',P3)';
   Kx=inv(T)*[zeros(na+nc-k,(m0+md)),Kx];
end
   
Bm=Bs(:,1:(m0+md));
if (m0+md)~=0 & na~=0,
   Bm(1:na,:)=Bm(1:na,:)+Kx(nc+1:nc+na,:);
end;
if (m0+md)~=0 & nc~=0,
   Bm(na+nb+1:na+nb+nc,:)=Bm(na+nb+1:na+nb+nc,:)+Kx(1:nc,:);
end;
Bm=Gs*Bm*invGm;
Dm=zeros(p,m0+md);
if m0~=0,
   Dm(1:m0,1:m0)=eye(m0);
end;
Dm=Go*Dm*invGm;

Av=Ax-Kx*Cx;
Bv=Bx-Kx*Dx;
Cv=inv(invGm)*Cx;
Dv=inv(invGm)*Dx;