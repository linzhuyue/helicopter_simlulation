function [Am,Bm,Cm,Dm,Av,Bv,Cv,Dv]=mpfact(A,B,C,D,tol)

%MPFACT  Minimum-Phase/All-Pass Factorization of Continuous-time Systems
%
%        [Am,Bm,Cm,Dm,Av,Bv,Cv,Dv] = mpfact(A,B,C,D)
%
%        calculates a minimum-phase/all-pass factorization for a
%        detectable system (A,B,C,D) with transfer function matrix
%        G(s), where [B' D'] and [C D] are assumed to be of full
%        rank, and (A,B,C,D) has no invariant zeros on the jw axis.
%
%        The minimum-phase/all-pass factorization is given as
%
%               G(s) = Gm(s) V(s)
%
%        where
%
%               Gm(s) = Cm (sI - Am)^{-1} Bm + Dm
%
%        is of minimum-phase and left invertible, and
%
%               V(s) = Cv (sI - Av)^{-1} Bv + Dv
%
%        is an all-pass factor satisfying V(s) V'(-s) = I.
%
%        See also IOFACT, GCFACT and DMPFACT.

if nargin==4,
   tol=1e-10;
end

pp=size(C,1);
mm=size(B,2);
if (mm==0)|(pp==0)
   disp(' ')
   disp('Warning: B or/and C is empty...');
   Am=[];Bm=[];Cm=[];Dm=[];Av=[];Bv=[];Cv=[];Dv=[];
   return;
end

m=rank([B;D],tol);
p=rank([C,D],tol);
if (mm~=m)|(pp~=p)
   disp(' ')
   disp('Warning: [B'' D''] or/and [C D] is not full rank...');
   Am=[];Bm=[];Cm=[];Dm=[];Av=[];Bv=[];Cv=[];Dv=[];
   return;
end

n=size(A,1);
[at,t1,t1,t1,ks]=obsvf(A,0*B,C,tol);
ks=sum(ks);
t1=eig(at(1:n-ks,1:n-ks));
if any(real(t1)>=0)
   disp('Warning: The system is not detectable...');
   Am=[];Bm=[];Cm=[];Dm=[];Av=[];Bv=[];Cv=[];Dv=[];
   return;
end

[As,Bs,Cs,Ds,Gs,Go,Gi,dims,lv,rv,qv,m0]=scb(A,B,C,D,tol);

nan=dims(1);na0=dims(2);nap=dims(3);
na=sum(dims(1:3));nb=dims(4);nc=dims(5);nd=dims(6);
md=length(qv);

if na0~=0
   disp('Warning: There are invariant zeros on the jw axis or unit circle (for discrete-time case).');
   Am=[];Bm=[];Cm=[];Dm=[];Av=[];Bv=[];Cv=[];Dv=[];
   return;
end;
[p,m]=size(C*B);

Am=A;Cm=C;
if (nap+nc)==0,
   Bm=B;Dm=D;Av=[];Bv=zeros(0,m);Cv=zeros(m,0);Dv=eye(m);
   return;
end

Ax=[As(na+nb+1:na+nb+nc,na+nb+1:na+nb+nc),As(na+nb+1:na+nb+nc,nan+1:nan+nap)];
Ax=[Ax;[zeros(nap,nc),As(nan+1:nan+nap,nan+1:nan+nap)]];
Bx=zeros(nap+nc,m);
if nc~=0,
   Bx(1:nc,m0+md+1:m)=Bs(na+nb+1:na+nb+nc,m0+md+1:m);
end
Bx=Bx*inv(Gi);
Bd=Bs(na+nb+nc+1:n,m0+1:m0+md);
Cx=[Cs(1:m0,na+nb+1:na+nb+nc),Cs(1:m0,nan+1:na)];
Cx=[Cx;Bd'*[As(na+nb+nc+1:n,na+nb+1:na+nb+nc),As(na+nb+nc+1:n,nan+1:na)]];
Dx=zeros(m0+md,m);
if (m0+md)~=0,
   Dx(:,1:(m0+md))=eye(m0+md);
end;
Dx=Dx*inv(Gi);
G=inv(Gi)*inv(Gi');
invGm=G(1:(m0+md),1:(m0+md));
if norm(invGm-eye((m0+md)))>tol*1e-4,
   invGm=sqrtm(invGm);
end;
AA=Ax-Bx*Dx'*inv(Dx*Dx')*Cx;
Q=Bx*Bx'-Bx*Dx'*inv(Dx*Dx')*Dx*Bx';
R=Cx'*inv(Dx*Dx')*Cx;
Px=are(AA',R,Q);
Kx=(Cx*Px+Dx*Bx')'*inv(Dx*Dx');
Bm=Bs(:,1:(m0+md));
if (m0+md)~=0 & nap~=0,
   Bm(nan+1:na,:)=Bm(nan+1:na,:)+Kx(nc+1:nc+nap,:);
end;
if (m0+md)~=0 & nc~=0,
   Bm(na+nb+1:na+nb+nc,:)=Bm(na+nb+1:na+nb+nc,:)+Kx(1:nc,:);
end;
Bm=Gs*Bm*(invGm);
Dm=zeros(p,(m0+md));
if m0~=0,
   Dm(1:m0,1:m0)=eye(m0);
end;
Dm=Go*Dm*(invGm);
Av=Ax-Kx*Cx;
Bv=Bx-Kx*Dx;
Cv=inv(invGm)*Cx;
Dv=inv(invGm)*Dx;