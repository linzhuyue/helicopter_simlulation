function [E,A,U,V,k]=zzschur2new(E,A,tol)

%Stability structural decomposition of the pair (E, A)
%
%     [E,A,U,V,k]=zzschur2new(EE,AA,tol)
%     gives the following transformation
%
%     U * E * V = EE = [E1  E*]         U * A * V = AA = [A1  A*]
%                      [ 0  E2],                         [ 0  A2],
%     
%     where eig(inv(E1)*A1) < -tol, and eig(inv(E2)*A2) >= -tol
%     size(E1,1)=k.

%EE, AA\in \R^{n\times n} with EE nonsingular

if nargin==2
   tol=1e-7;
end

n=size(E,1);
U(1:n,1:n)=eye(n);
V(1:n,1:n)=eye(n);
xnn=n;
XXE(1:xnn,1:xnn)=E(1:n,1:n);
XXA(1:xnn,1:xnn)=A(1:n,1:n);
k=0;
[VV(1:xnn,1:xnn),DD(1:xnn,1:xnn)]=eig(XXA(1:xnn,1:xnn),XXE(1:xnn,1:xnn));

while (xnn>0)
   xxkk=0;
   for i=1:xnn
      if (real(DD(i,i))<-tol)
         xxkk=1;
         i=xnn;
      end
   end
   if (xxkk==0)
      xnn=0;
   end
   if (xxkk==1)
      xkkk=0;
      [XVV(1:xnn,1:xnn),XDD(1:xnn,1:xnn)]=eig(XXA(1:xnn,1:xnn),XXE(1:xnn,1:xnn));
      xkk=0;
      i=1;
      while i<=xnn, 
         if real(XDD(i,i))<-tol
            xkk=i;
            i=xnn;
         end
         i=i+1;
      end
      if xkk>0
         for i=(xkk+1):xnn
            if (real(XDD(i,i))<-tol)&(abs(real(XDD(i,i)))+tol<abs(real(XDD(xkk,xkk))))
               xkk=i;
            end
         end
         XQQ(1:xnn,1)=real(XVV(1:xnn,xkk));
         XQQ(1:xnn,2)=imag(XVV(1:xnn,xkk));
      end
      XXV(1:xnn,1:xnn)=eye(xnn);
      XXU(1:xnn,1:xnn)=eye(xnn);
      XAA(1:xnn,1:xnn)=XXA(1:xnn,1:xnn);
      XEE(1:xnn,1:xnn)=XXE(1:xnn,1:xnn);
      if xkk>0
         [XXV(1:xnn,1:xnn),XDD1(1:xnn,1:2),XP1(1:2,1:2)]=svd(XQQ(1:xnn,1:2));
         xkkk=rank(XDD1(1:xnn,1:2),tol);
         [XXU(1:xnn,1:xnn),XDD2(1:xnn,1:xkkk),XP2(1:xkkk,1:xkkk)]=svd(XXE(1:xnn,1:xnn)*XXV(1:xnn,1:xkkk));
         XXU(1:xnn,1:xnn)=XXU(1:xnn,1:xnn)';
         XAA(1:xnn,1:xnn)=XXU(1:xnn,1:xnn)*XXA(1:xnn,1:xnn)*XXV(1:xnn,1:xnn);
         XEE(1:xnn,1:xnn)=XXU(1:xnn,1:xnn)*XXE(1:xnn,1:xnn)*XXV(1:xnn,1:xnn);
         XAA((xkkk+1):xnn,1:xkkk)=0;
         XEE((xkkk+1):xnn,1:xkkk)=0;
      end
      k=k+xkkk;
      E((n-xnn+1):n,(n-xnn+1):n)=XEE(1:xnn,1:xnn);
      A((n-xnn+1):n,(n-xnn+1):n)=XAA(1:xnn,1:xnn);
      U((n-xnn+1):n,1:n)=XXU(1:xnn,1:xnn)*U((n-xnn+1):n,1:n);
      V(1:n,(n-xnn+1):n)=V(1:n,(n-xnn+1):n)*XXV(1:xnn,1:xnn);
      E(1:(n-xnn),(n-xnn+1):n)=E(1:(n-xnn),(n-xnn+1):n)*XXV(1:xnn,1:xnn);
      A(1:(n-xnn),(n-xnn+1):n)=A(1:(n-xnn),(n-xnn+1):n)*XXV(1:xnn,1:xnn);
      XXE(1:(xnn-xkkk),1:(xnn-xkkk))=XEE((xkkk+1):xnn,(xkkk+1):xnn);
      XXA(1:(xnn-xkkk),1:(xnn-xkkk))=XAA((xkkk+1):xnn,(xkkk+1):xnn);
      xnn=xnn-xkkk;
      if xnn>0
         [VV(1:xnn,1:xnn),DD(1:xnn,1:xnn)]=eig(XXA(1:xnn,1:xnn),XXE(1:xnn,1:xnn));
      end
   end
end