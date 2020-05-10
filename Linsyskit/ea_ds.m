function [PEQ,PAQ,P,Q,n1,n2]=ea_ds(E,A,tol)

%EA_DS  Decomposition for a Pair of Square Matrices (E,A)
%
%        [Et,At,P,Q,n1,n2] = ea_ds(E,A)
%
%        decomposes (E,A), related to descriptor systems, i.e.,
%                               .
%                             E x = A x + ...
%
%        into the following special form:
%
%                     Et = P E Q = [ I   0 ]  n1
%                                  [ 0   N ]  n2
%
%                     At = P A Q = [ A1  0 ]  n1
%                                  [ 0   I ]  n2
%
%        where N is a nilpotent matrix, P and Q are nonsingular.
%
%        See also SD_DS.

if nargin==2
   tol=1e-8;
end

n=size(A,1);

if isempty(A)
   PEQ=[];PAQ=[];P=[];Q=[];n1=[];n2=[];
   return
end

[EE,AA,P,Q,nn,ll]=zzguptri2(E,A,tol);
n1=nn(2);n2=nn(3);

if (n1+n2)~=n
   disp('  Check det(s*E-A), might be it always equal to zero........')
   PEQ=E;PAQ=A;P=eye(n);Q=eye(n);n1=[];n2=[];
   return
end

if n1==0
   [u,s,v]=svd(A); tt=[];
   for kk=1:n
      tt(kk,kk)=sqrt(s(kk,kk));
   end
   P=u*tt;
   Q=tt*v;
   PEQ=P*E*Q;
   PAQ=eye(n);
   return
end

if n2==0
   [u,s,v]=svd(E); tt=[];
   for kk=1:n
      tt(kk,kk)=sqrt(s(kk,kk));
   end
   P=u*tt;
   Q=tt*v;
   PEQ=eye(n);
   PAQ=P*A*Q;
   return
end

E11=EE(1:n1,1:n1);E12=EE(1:n1,n1+1:n);E22=EE(n1+1:n,n1+1:n);
A11=AA(1:n1,1:n1);A12=AA(1:n1,n1+1:n);A22=AA(n1+1:n,n1+1:n);

[Y,X]=zzgensyl(E11,E22,-E12,A11,A22,-A12);
Pt=eye(n);Pt(1:n1,n1+1:n)=X;
Qt=eye(n);Qt(1:n1,n1+1:n)=Y;

P=Pt*P;
Q=Q*Qt;
te=blkdiag(inv(E11),eye(n2));
ta=blkdiag(eye(n1),inv(A22));

t(1)=max(cond(te*ta*P),cond(Q));
t(2)=max(cond(te*P),cond(Q*ta));
t(3)=max(cond(ta*P),cond(Q*te));
t(4)=max(cond(P),cond(Q*ta*te));

[y,h]=sort(t);
switch h(1)
   case 1
      P=te*ta*P;
      A11=inv(E11)*A11;
      E22=inv(A22)*E22;
   case 2
      P=te*P;Q=Q*ta;
      A11=inv(E11)*A11;
      E22=E22*inv(A22);
   case 3
      P=ta*P;Q=Q*te;
      A11=A11*inv(E11);
      E22=inv(A22)*E22;
   case 4
      Q=Q*ta*te;
      A11=A11*inv(E11);
      E22=E22*inv(A22);
end

PEQ=blkdiag(eye(n1),E22);
PAQ=blkdiag(A11,eye(n2));
%cond(P),cond(Q)

err_E=norm(PEQ-P*E*Q);err_A=norm(PAQ-P*A*Q);