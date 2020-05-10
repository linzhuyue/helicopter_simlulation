function [U,V,E,A,xn12,xl12]=guptri1(E,A,xn,xl,tol)

%GUPTRI1  [U,V,EE,EA,xn12,xl12]=zzguptri1(E,A,xn,xl,tol)
%
%     gives the following transformation:
%
%                                           xl12     xl-xl12
%     U*(sE(1:xn,1:xl)-A(1:xn,1:xl))*V = [sE11-A11  sE12-A12]   xn12
%                                        [    0     sE22-A22]   xn-xn12
%
%     where E11 is of full-row rank;
%           sE22-A22 is of full-column rank;

%U(1:xn,1:xn)=eye(xn);
%V(1:xl,1:xl)=eye(xl);
U=eye(xn);
V=eye(xl);
xn12=xn;
xl12=xl;
t=rank(E(1:xn12,1:xl12),tol);
while (t<xn12)
   [Q(1:xn12,1:xn12),R(1:xn12,1:xl12),P(1:xl12,1:xl12)]=svd(E(1:xn12,1:xl12));
   Q(1:xn12,1:xn12)=Q(1:xn12,1:xn12)';
   P(1:xl12,1:xl12)=P(1:xl12,1:xl12)';
   U(1:xn12,1:xn)=Q(1:xn12,1:xn12)*U(1:xn12,1:xn);
   R((t+1):xn12,1:xl12)=0;
   E(1:t,1:xl12)=R(1:t,1:xl12)*P(1:xl12,1:xl12);
   E((t+1):xn12,1:xl12)=0;
   E(1:xn12,(xl12+1):xl)=Q(1:xn12,1:xn12)*E(1:xn12,(xl12+1):xl);
   A(1:xn12,1:xl)=Q(1:xn12,1:xn12)*A(1:xn12,1:xl);
   [Q(1:(xn12-t),1:(xn12-t)),R(1:(xn12-t),1:xl12),P(1:xl12,1:xl12)]=svd(A((t+1):xn12,1:xl12));
   
   s=rank(R(1:(xn12-t),1:xl12),tol);
   R((s+1):(xn12-t),1:xl12)=0;
   XII(1:xl12,1:xl12)=0;
   XII(1:s,(xl12-s+1):xl12)=eye(s);
   XII((s+1):xl12,1:(xl12-s))=eye(xl12-s);
   P(1:xl12,1:xl12)=P(1:xl12,1:xl12)*XII(1:xl12,1:xl12);
   A((t+1):xn12,1:xl12)=Q(1:(xn12-t),1:(xn12-t))*R(1:(xn12-t),1:xl12)*XII(1:xl12,1:xl12);
   A((t+1):xn12,1:(xl12-s))=0;
   E(1:t,1:xl12)=E(1:t,1:xl12)*P(1:xl12,1:xl12);
   A(1:t,1:xl12)=A(1:t,1:xl12)*P(1:xl12,1:xl12);
   V(1:xl,1:xl12)=V(1:xl,1:xl12)*P(1:xl12,1:xl12);
   xn12=t;
   xl12=xl12-s;
   t=rank(E(1:xn12,1:xl12),tol);
end