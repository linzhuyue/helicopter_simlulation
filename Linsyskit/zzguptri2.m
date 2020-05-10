function [E,A,U,V,nn,ll]=zzguptri2(E,A,tol)

%     [EE,AA,U,V,nn,ll,na0]=zzguptri2(E,A,tol);
%
%                     ll(1)     ll(2)     ll(3)     ll(4)
%                  [sE11-A11  sE12-A12  sE13-A13  sE14-A14]   nn(1)
%     U*(sE-A)*V = [    0     sE22-A22  sE23-A23  sE24-A24]   nn(2)
%                  [    0         0     sE33-A33  sE34-A34]   nn(3)
%                  [    0         0         0     sE44-A44]   nn(4)
%
%     where E11 is of full-row rank;
%           sE11-A11 is of full-row rank;
%           E22 is of nonsigular;
%           sE33-A33 is of nonsigular;
%           E44 is of full-column rank;
%           sE44-A44 is of full-column rank.

[nn,ll]=size(E);
XE=E;
XA=A;

xn=nn;
xl=ll;
[XU,XV,XE,XA,xn12,xl12]=zzguptri1(XE,XA,xn,xl,tol);

U(1:xn,1:xn)=XU(1:xn,1:xn);
V(1:xl,1:xl)=XV(1:xl,1:xl);
E(1:xn,1:xl)=XE(1:xn,1:xl);
A(1:xn,1:xl)=XA(1:xn,1:xl);

ttss=xn12;
sstt=xl12;

if (xn12*xl12>0)
   xn=xl12;
   xl=xn12;
   XE(1:xn,1:xl)=E(1:xn12,1:xl12)';
   XA(1:xn,1:xl)=A(1:xn12,1:xl12)';
   
   [XU,XV,XE,XA,xn12,xl12]=zzguptri1(XE,XA,xn,xl,tol);
   
   XZ(1:xl,1:xl)=0;
   XZ((xl-xl12+1):xl,1:xl12)=eye(xl12);
   XZ(1:(xl-xl12),(xl12+1):xl)=eye(xl-xl12);
   XW(1:xn,1:xn)=0;
   XW(1:xn12,(xn-xn12+1):xn)=eye(xn12);
   XW((xn12+1):xn,1:(xn-xn12))=eye(xn-xn12);
   E(1:xl,1:xn)=XZ(1:xl,1:xl)*XE(1:xn,1:xl)'*XW(1:xn,1:xn);
   A(1:xl,1:xn)=XZ(1:xl,1:xl)*XA(1:xn,1:xl)'*XW(1:xn,1:xn);
   nn1=xl-xl12;
   nn2=xl12;
   ll1=xn-xn12;
   ll2=xn12;
   U(1:xl,1:nn)=XZ(1:xl,1:xl)*(XV(1:xl,1:xl)')*U(1:xl,1:nn);
   V(1:ll,1:xn)=V(1:ll,1:xn)*(XU(1:xn,1:xn)')*XW(1:xn,1:xn);
   E(1:xl,(xn+1):ll)=XZ(1:xl,1:xl)*(XV(1:xl,1:xl)')*E(1:xl,(xn+1):ll);
   A(1:xl,(xn+1):ll)=XZ(1:xl,1:xl)*(XV(1:xl,1:xl)')*A(1:xl,(xn+1):ll);
end

if (ttss==0)
   nn1=0;
   nn2=0;
   ll1=sstt;
   ll2=0;
end

if (sstt==0)
   ll1=0;
   ll2=0;
   nn1=0;
   nn2=0;
   ttss=0;
end

xn=nn-ttss;
xl=ll-sstt;

if (xn*xl>0)
   XE(1:xn,1:xl)=A((nn-xn+1):nn,(ll-xl+1):ll);
   XA(1:xn,1:xl)=E((nn-xn+1):nn,(ll-xl+1):ll);

   [XU,XV,XE,XA,xn12,xl12]=zzguptri1(XE,XA,xn,xl,tol);

   nn3=xn12;
   nn4=xn-xn12;
   ll3=xl12;
   ll4=xl-ll3;
   E((nn-xn+1):nn,(ll-xl+1):ll)=XA(1:xn,1:xl);
   A((nn-xn+1):nn,(ll-xl+1):ll)=XE(1:xn,1:xl);
   U((nn-xn+1):nn,1:nn)=XU(1:xn,1:xn)*U((nn-xn+1):nn,1:nn);
   V(1:ll,(ll-xl+1):ll)=V(1:ll,(ll-xl+1):ll)*XV(1:xl,1:xl);
   E(1:(nn-xn),(ll-xl+1):ll)=E(1:(nn-xn),(ll-xl+1):ll)*XV(1:xl,1:xl);
   A(1:(nn-xn),(ll-xl+1):ll)=A(1:(nn-xn),(ll-xl+1):ll)*XV(1:xl,1:xl);
end

if (xn==0)
   nn3=0;
   nn4=0;
   ll3=0;
   ll4=0;
end

if (xl==0)
   ll3=0;
   ll4=0;
   nn3=0;
   nn4=xn;
end

nn=[nn1 nn2 nn3 nn4];
ll=[ll1 ll2 ll3 ll4];