function [U56,V56,QB56,CP56,A56,B56,C56,Ind]=zzinfzerogsvd(E56,A56,B56,C56,n5,n6,tol)

%Infzerogsvd [U,V,Gi,Go,AA,BB,CC,Ind]=zzinfzerogsvd(E,A,B,C,n5,n6,tol)
%    
%      Extraction of Infinite Zero Strucuture.
%
%      Input:
%             System (E,A,B,C) has only infinite zeros.
%             n5 is the rank of C and B;
%             n6=size(A,1)-n5;
%             B(1:n6,:)=0; C(:,1:n6)=0;
%
%      Output: 
%             s*I-AA = U*(s*E-A)*V
%                 BB = U*B*Gi
%                 CC = Go*C*V
%           where (AA,BB,CC) are in the form of (Ad,Bd,Cd).
%           The orders of infinite zeros are descending.

nn=n5+n6;
Ind=[];
p=[];
QB(1:nn,1:nn)=eye(nn);
CP(1:nn,1:nn)=eye(nn);

if n5>0
   k=0;
   t=0;   
   E5(1:n5,1:n5)=E56(1:n5,1:n5);
   A5(1:n5,1:n5)=A56(1:n5,1:n5);
   U56(1:nn,1:nn)=eye(nn);
   V56(1:nn,1:nn)=eye(nn);
   U5(1:n5,1:n5)=eye(n5);
   V5(1:n5,1:n5)=eye(n5);
   while t<n5
      k=k+1;
      tt=n5-t;
      [UX(1:tt,1:tt),VX(1:tt,1:tt),XX(1:tt,1:tt),S5(1:tt,1:tt),SX(1:tt,1:tt)]=gsvd(E5((t+1):n5,(t+1):n5),A5((t+1):n5,(t+1):n5));
      r=rank(S5(1:tt,1:tt),tol);
      ii=tt-r;
      if ii==0
         ii=1;
         for jk=2:tt
            if S5(jk,jk)<(S5(1,1)+10.0*tol)
               ii=jk;
            end
         end 
      end
      U5((t+1):n5,1:n5)=VX(1:tt,1:tt)'*U5((t+1):n5,1:n5);
      V5(1:n5,(t+1):n5)=V5(1:n5,(t+1):n5)*inv(XX(1:tt,1:tt)');
      A5((t+1):n5,(t+1):n5)=SX(1:tt,1:tt);
      E5((t+1):n5,(t+1):n5)=VX(1:tt,1:tt)'*UX(1:tt,1:tt)*S5(1:tt,1:tt);
      E5(1:t,(t+1):n5)=E5(1:t,(t+1):n5)*inv(XX(1:tt,1:tt)');
      E5((t+1):n5,(t+1):(t+ii))=0;
      t=t+ii;
      p(k)=ii;
   end
   t2=n5;
   for i=1:(k-1)
      t1=t2-p(k+1-i);
      t3=p(k-i);
      t4=p(k+1-i);
      [Q(1:t3,1:t3),R(1:t3,1:t4),EX(1:t4,1:t4)]=qr(E5((t1-t3+1):t1,(t1+1):t2));
      Q(1:t3,1:t3)=Q(1:t3,1:t3)';
      EX(1:t4,1:t4)=EX(1:t4,1:t4)';
      U5((t1-t3+1):t1,1:n5)=Q(1:t3,1:t3)*U5((t1-t3+1):t1,1:n5);
      A5((t1-t3+1):t1,(t1-t3+1):t1)=Q(1:t3,1:t3)*A5((t1-t3+1):t1,(t1-t3+1):t1);
      E5((t1-t3+1):t1,(t2+1):n5)=Q(1:t3,1:t3)*E5((t1-t3+1):t1,(t2+1):n5);
      E5((t1-t3+t4+1):t1,(t1+1):t2)=0;
      E5((t1-t3+1):(t1-t3+t4),(t1+1):t2)=R(1:t4,1:t4)*EX(1:t4,1:t4);
      t2=t1;
   end
   V5(1:n5,1:n5)=V5(1:n5,1:n5)*inv(A5(1:n5,1:n5));
   E5(1:n5,1:n5)=E5(1:n5,1:n5)*inv(A5(1:n5,1:n5));
   A5(1:n5,1:n5)=eye(n5);
   EA(1:n5,1:n5)=E5(1:n5,1:n5);
   AA(1:n5,1:n5)=A5(1:n5,1:n5);
   UA(1:n5,1:n5)=U5(1:n5,1:n5);
   VA(1:n5,1:n5)=V5(1:n5,1:n5);
   U5(1:n5,1:n5)=eye(n5);
   V5(1:n5,1:n5)=eye(n5);
   t2=n5;
   for ih=2:(k-1)
      i=k+1-ih;
      t4=p(i+1);
      t1=t2-t4;
      t3=t1-p(i);
      E5(1:t3,(t3+1):(t3+t4))=E5(1:t3,(t3+1):(t3+t4))*E5((t3+1):(t3+t4),(t1+1):t2);
      V5(1:n5,(t3+1):(t3+t4))=V5(1:n5,(t3+1):(t3+t4))*E5((t3+1):(t3+t4),(t1+1):t2);
      U5((t3+1):(t3+t4),1:n5)=inv(E5((t3+1):(t3+t4),(t1+1):t2))*U5((t3+1):(t3+t4),1:n5);
      E5((t3+1):(t3+t4),(t1+1):t2)=eye(t4);
      Q(1:t3,1:t4)=-E5(1:t3,(t1+1):t2);
      E5(1:t3,(t1+1):t2)=0;
      U5(1:t3,1:n5)=U5(1:t3,1:n5)+Q(1:t3,1:t4)*U5((t3+1):(t3+t4),1:n5);
      E5(1:t3,(t3+1):(t3+t4))=E5(1:t3,(t3+1):(t3+t4))-E5(1:t3,1:t3)*Q(1:t3,1:t4);
      V5(1:n5,(t3+1):(t3+t4))=V5(1:n5,(t3+1):(t3+t4))-V5(1:n5,1:t3)*Q(1:t3,1:t4);
      t2=t1;
   end
   if k>1
      U5(1:p(2),1:n5)=inv(E5(1:p(2),(p(1)+1):(p(1)+p(2))))*U5(1:p(2),1:n5);
      V5(1:n5,1:p(2))=V5(1:n5,1:p(2))*E5(1:p(2),(p(1)+1):(p(1)+p(2)));
      E5(1:p(2),(p(1)+1):(p(1)+p(2)))=eye(p(2));
   end
   U56(1:n5,1:nn)=U5(1:n5,1:n5)*UA(1:n5,1:n5)*U56(1:n5,1:nn);
   V56(1:nn,1:n5)=V56(1:nn,1:n5)*VA(1:n5,1:n5)*V5(1:n5,1:n5);
   E56(1:n5,1:n5)=E5(1:n5,1:n5);
   A56(1:n5,1:n5)=A5(1:n5,1:n5);
   E56(1:n5,(n5+1):nn)=U5(1:n5,1:n5)*UA(1:n5,1:n5)*E56(1:n5,(n5+1):nn);
   A56(1:n5,(n5+1):nn)=U5(1:n5,1:n5)*UA(1:n5,1:n5)*A56(1:n5,(n5+1):nn);
   E56((n5+1):nn,1:n5)=E56((n5+1):nn,1:n5)*VA(1:n5,1:n5)*V5(1:n5,1:n5);
   A56((n5+1):nn,1:n5)=A56((n5+1):nn,1:n5)*VA(1:n5,1:n5)*V5(1:n5,1:n5);
   t=0;
   tt=p(1);
   for i=2:k
      A56((n5+1):nn,1:nn)=A56((n5+1):nn,1:nn)-E56((n5+1):nn,(tt+1):(tt+p(i)))*A56((t+1):(t+p(i)),1:nn);
      E56((n5+1):nn,(n5+1):nn)=E56((n5+1):nn,(n5+1):nn)-E56((n5+1):nn,(tt+1):(tt+p(i)))*E56((t+1):(t+p(i)),(n5+1):nn);
      U56((n5+1):nn,1:nn)=U56((n5+1):nn,1:nn)-E56((n5+1):nn,(tt+1):(tt+p(i)))*U56((t+1):(t+p(i)),1:nn);
      t=t+p(i-1);
      tt=tt+p(i);
   end
   E56((n5+1):nn,(p(1)+1):n5)=0;
   [Q(1:n6,1:n6),R(1:n6,1:p(1)),EX(1:p(1),1:p(1))]=qr(E56((n5+1):nn,1:p(1)));
   Q(1:n6,1:n6)=Q(1:n6,1:n6)';
   U56((n5+1):nn,1:nn)=Q(1:n6,1:n6)*U56((n5+1):nn,1:nn);
   A56((n5+1):nn,1:nn)=Q(1:n6,1:n6)*A56((n5+1):nn,1:nn);
   E56((n5+1):nn,1:nn)=Q(1:n6,1:n6)*E56((n5+1):nn,1:nn);
   B56((n5+1):nn,1:n6)=Q(1:n6,1:n6)*B56((n5+1):nn,1:n6);
   U56((n5+1):(n5+p(1)),1:nn)=EX(1:p(1),1:p(1))*inv(R(1:p(1),1:p(1)))*U56((n5+1):(n5+p(1)),1:nn);
   A56((n5+1):(n5+p(1)),1:nn)=EX(1:p(1),1:p(1))*inv(R(1:p(1),1:p(1)))*A56((n5+1):(n5+p(1)),1:nn);
   E56((n5+1):(n5+p(1)),(p(1)+1):nn)=EX(1:p(1),1:p(1))*inv(R(1:p(1),1:p(1)))*E56((n5+1):(n5+p(1)),(p(1)+1):nn);
   B56((n5+1):(n5+p(1)),1:n6)=EX(1:p(1),1:p(1))*inv(R(1:p(1),1:p(1)))*B56((n5+1):(n5+p(1)),1:n6);
   E56((n5+1):nn,1:p(1))=0;
   E56((n5+1):(n5+p(1)),1:p(1))=eye(p(1));
   V56(1:nn,(n5+1):nn)=V56(1:nn,(n5+1):nn)-V56(1:nn,1:p(1))*E56((n5+1):(n5+p(1)),(n5+1):nn);
   A56(1:nn,(n5+1):nn)=A56(1:nn,(n5+1):nn)-A56(1:nn,1:p(1))*E56((n5+1):(n5+p(1)),(n5+1):nn);
   E56((n5+1):(n5+p(1)),(n5+1):nn)=0;
   t=0;
   tt=p(1);
   for i=2:k
      V56(1:nn,(n5+1):nn)=V56(1:nn,(n5+1):nn)-V56(1:nn,(tt+1):(tt+p(i)))*E56((t+1):(t+p(i)),(n5+1):nn);
      A56(1:nn,(n5+1):nn)=A56(1:nn,(n5+1):nn)-A56(1:nn,(tt+1):(tt+p(i)))*E56((t+1):(t+p(i)),(n5+1):nn);
      E56((t+1):(t+p(i)),(n5+1):nn)=0;
      t=t+p(i-1);
      tt=tt+p(i);
   end
   X(1:p(1),1:n6)=0;
   t=0;
   tt=0;
   for i=1:(k-1)
      X((t+1):(t+p(i)-p(i+1)),1:n6)=E56((tt+p(i+1)+1):(tt+p(i)),(n5+1):nn);
      t=t+p(i)-p(i+1);
      tt=tt+p(i);
   end
   X((t+1):p(1),1:n6)=E56((tt+1):n5,(n5+1):nn);
   [Q(1:n6,1:n6),R(1:n6,1:p(1)),EX(1:p(1),1:p(1))]=qr(X(1:p(1),1:n6)');
   V56(1:nn,(n5+1):nn)=V56(1:nn,(n5+1):nn)*Q(1:n6,1:n6);
   A56(1:nn,(n5+1):nn)=A56(1:nn,(n5+1):nn)*Q(1:n6,1:n6);
   E56((n5+p(1)+1):nn,(n5+1):nn)=E56((n5+p(1)+1):nn,(n5+1):nn)*Q(1:n6,1:n6);
   C56(1:n6,(n5+1):nn)=C56(1:n6,(n5+1):nn)*Q(1:n6,1:n6);
   V56(1:nn,(n5+1):(n5+p(1)))=V56(1:nn,(n5+1):(n5+p(1)))*inv(R(1:p(1),1:p(1))')*EX(1:p(1),1:p(1))';
   A56(1:nn,(n5+1):(n5+p(1)))=A56(1:nn,(n5+1):(n5+p(1)))*inv(R(1:p(1),1:p(1))')*EX(1:p(1),1:p(1))';
   E56((n5+p(1)+1):nn,(n5+1):(n5+p(1)))=E56((n5+p(1)+1):nn,(n5+1):(n5+p(1)))*inv(R(1:p(1),1:p(1))')*EX(1:p(1),1:p(1))';
   C56(1:n6,(n5+1):(n5+p(1)))=C56(1:n6,(n5+1):(n5+p(1)))*inv(R(1:p(1),1:p(1))')*EX(1:p(1),1:p(1))';
   X(1:p(1),1:n6)=0;
   X(1:p(1),1:p(1))=eye(p(1));
   t=n5;
   tt=0;
   for i=1:(k-1)
      E56((tt+p(i+1)+1):(tt+p(i)),(n5+1):nn)=0;
      E56((tt+p(i+1)+1):(tt+p(i)),(t+1):(t+p(i)-p(i+1)))=eye(p(i)-p(i+1));
      t=t+p(i)-p(i+1);
      tt=tt+p(i);
   end   
   E56((tt+1):n5,(n5+1):nn)=0;
   E56((tt+1):n5,(t+1):(n5+p(1)))=eye(p(k));
   U56((n5+p(1)+1):nn,1:nn)=inv(E56((n5+p(1)+1):nn,(n5+p(1)+1):nn))*U56((n5+p(1)+1):nn,1:nn);
   A56((n5+p(1)+1):nn,1:nn)=inv(E56((n5+p(1)+1):nn,(n5+p(1)+1):nn))*A56((n5+p(1)+1):nn,1:nn);
   B56((n5+p(1)+1):nn,1:n6)=inv(E56((n5+p(1)+1):nn,(n5+p(1)+1):nn))*B56((n5+p(1)+1):nn,1:n6);
   E56((n5+p(1)+1):nn,(n5+1):(n5+p(1)))=inv(E56((n5+p(1)+1):nn,(n5+p(1)+1):nn))*E56((n5+p(1)+1):nn,(n5+1):(n5+p(1)));
   E56((n5+p(1)+1):nn,(n5+p(1)+1):nn)=eye(n6-p(1));
   V56(1:nn,(n5+1):(n5+p(1)))=V56(1:nn,(n5+1):(n5+p(1)))-V56(1:nn,(n5+p(1)+1):nn)*E56((n5+p(1)+1):nn,(n5+1):(n5+p(1)));
   A56(1:nn,(n5+1):(n5+p(1)))=A56(1:nn,(n5+1):(n5+p(1)))-A56(1:nn,(n5+p(1)+1):nn)*E56((n5+p(1)+1):nn,(n5+1):(n5+p(1)));
   C56(1:n6,(n5+1):(n5+p(1)))=C56(1:n6,(n5+1):(n5+p(1)))-C56(1:n6,(n5+p(1)+1):nn)*E56((n5+p(1)+1):nn,(n5+1):(n5+p(1)));
   E56((n5+p(1)+1):nn,(n5+1):(n5+p(1)))=0;
   t=0;
   for i=1:k
      ik=(k+1)-i;
      IIND(1,i)=ik;
      IIND(2,i)=p(ik)-t;
      t=p(ik);
   end
   J=0;
   for i=1:k
      for j=1:IIND(2,i)
         J=J+1;
         Ind(J)=IIND(1,i);
      end
   end
   J=1;
   ex(1:n5,1:n5)=eye(n5);
   for i=1:p(1)
      VV(1:n5,J:J)=ex(1:n5,i:i);
      t=0;
      for j=1:(Ind(i)-1)
          J=J+1;
          t=t+p(j);
          VV(1:n5,J:J)=ex(1:n5,(i+t):(i+t));
      end
      J=J+1;
   end
   V5(1:n5,1:n5)=V5(1:n5,1:n5)*VV(1:n5,1:n5);
   U5(1:n5,1:n5)=VV(1:n5,1:n5)'*U5(1:n5,1:n5);
   E5(1:n5,1:n5)=VV(1:n5,1:n5)'*E5(1:n5,1:n5)*VV(1:n5,1:n5);
   V56(1:nn,1:n5)=V56(1:nn,1:n5)*VV(1:n5,1:n5);
   E56(1:n5,1:n5)=VV(1:n5,1:n5)'*E56(1:n5,1:n5)*VV(1:n5,1:n5);
   E56((n5+1):nn,1:n5)=E56((n5+1):nn,1:n5)*VV(1:n5,1:n5);
   A56((n5+1):nn,1:n5)=A56((n5+1):nn,1:n5)*VV(1:n5,1:n5);
   E56(1:n5,(n5+1):nn)=VV(1:n5,1:n5)'*E56(1:n5,(n5+1):nn);
   A56(1:n5,(n5+1):nn)=VV(1:n5,1:n5)'*A56(1:n5,(n5+1):nn);
   U56(1:n5,1:nn)=VV(1:n5,1:n5)'*U56(1:n5,1:nn);
   VC(1:n5,1:n5)=0;
   [i,Jj]=size(Ind);% Jj is the number of the jordan block in sE5-A5.
   t=0;
   for i=1:Jj
      for j=1:Ind(i)
         VC(t+j,(t+Ind(i)+1-j))=1;
      end
      t=t+Ind(i);
   end
   V5(1:n5,1:n5)=V5(1:n5,1:n5)*VC(1:n5,1:n5);
   U5(1:n5,1:n5)=VC(1:n5,1:n5)*U5(1:n5,1:n5);
   E5(1:n5,1:n5)=VC(1:n5,1:n5)*E5(1:n5,1:n5)*VC(1:n5,1:n5);
   V56(1:nn,1:n5)=V56(1:nn,1:n5)*VC(1:n5,1:n5);
   E56(1:n5,1:n5)=VC(1:n5,1:n5)*E56(1:n5,1:n5)*VC(1:n5,1:n5);
   E56((n5+1):nn,1:n5)=E56((n5+1):nn,1:n5)*VC(1:n5,1:n5);
   A56((n5+1):nn,1:n5)=A56((n5+1):nn,1:n5)*VC(1:n5,1:n5);
   E56(1:n5,(n5+1):nn)=VC(1:n5,1:n5)*E56(1:n5,(n5+1):nn);
   A56(1:n5,(n5+1):nn)=VC(1:n5,1:n5)*A56(1:n5,(n5+1):nn);
   U56(1:n5,1:nn)=VC(1:n5,1:n5)*U56(1:n5,1:nn);
   X(1:p(1),1:p(1))=0;
   Y(1:p(1),1:p(1))=0;
   tt=0;
   for i=1:Jj
      X(i,1:p(1))=E56((tt+1),(n5+1):(n5+p(1)));
      Y(1:p(1),i)=E56((n5+1):(n5+p(1)),(tt+Ind(i)));
      tt=tt+Ind(i);
   end
   E56(1:n5,(n5+1):(n5+p(1)))=E56(1:n5,(n5+1):(n5+p(1)))*inv(X(1:p(1),1:p(1)));
   A56(1:nn,(n5+1):(n5+p(1)))=A56(1:nn,(n5+1):(n5+p(1)))*inv(X(1:p(1),1:p(1)));
   V56(1:nn,(n5+1):(n5+p(1)))=V56(1:nn,(n5+1):(n5+p(1)))*inv(X(1:p(1),1:p(1)));
   C56(1:n6,(n5+1):(n5+p(1)))=C56(1:n6,(n5+1):(n5+p(1)))*inv(X(1:p(1),1:p(1)));
   E56((n5+1):(n5+p(1)),1:n5)=inv(Y(1:p(1),1:p(1)))*E56((n5+1):(n5+p(1)),1:n5);
   A56((n5+1):(n5+p(1)),1:nn)=inv(Y(1:p(1),1:p(1)))*A56((n5+1):(n5+p(1)),1:nn);
   U56((n5+1):(n5+p(1)),1:nn)=inv(Y(1:p(1),1:p(1)))*U56((n5+1):(n5+p(1)),1:nn);
   B56((n5+1):(n5+p(1)),1:n6)=inv(Y(1:p(1),1:p(1)))*B56((n5+1):(n5+p(1)),1:n6);
   QB56(1:n6,1:n6)=inv(B56((n5+1):nn,1:n6));
   CP56(1:n6,1:n6)=inv(C56(1:n6,(n5+1):nn));
   B56((n5+1):nn,1:n6)=eye(n6);
   C56(1:n6,(n5+1):nn)=eye(n6);
   t=0;
   for i=1:Jj
      X(1:nn,1:(n5+i-1-t))=E56(1:nn,(t+1):(n5+i-1));
      Y(1:nn,1:(n5+i-1-t))=A56(1:nn,(t+1):(n5+i-1));
      Z(1:nn,1:(n5+i-1-t))=V56(1:nn,(t+1):(n5+i-1));
      W(1:n6,1:(n5+i-1-t))=C56(1:n6,(t+1):(n5+i-1));
      E56(1:nn,t+1)=E56(1:nn,n5+i);
      E56(1:nn,(t+2):(n5+i))=X(1:nn,1:(n5+i-1-t));
      A56(1:nn,t+1)=A56(1:nn,n5+i);
      A56(1:nn,(t+2):(n5+i))=Y(1:nn,1:(n5+i-1-t));
      V56(1:nn,t+1)=V56(1:nn,n5+i);
      V56(1:nn,(t+2):(n5+i))=Z(1:nn,1:(n5+i-1-t));
      C56(1:n6,t+1)=C56(1:n6,n5+i);
      C56(1:n6,(t+2):(n5+i))=W(1:n6,1:(n5+i-1-t));
      t=t+Ind(i)+1;
   end
   t=0;
   for i=1:Jj
      t=t+Ind(i);
      X(1:(n5+i-1-t),1:nn)=E56((t+1):(n5+i-1),1:nn);
      Y(1:(n5+i-1-t),1:nn)=A56((t+1):(n5+i-1),1:nn);
      Z(1:(n5+i-1-t),1:nn)=U56((t+1):(n5+i-1),1:nn);
      W(1:(n5+i-1-t),1:n6)=B56((t+1):(n5+i-1),1:n6);
      E56(t+1,1:nn)=E56(n5+i,1:nn);
      E56((t+2):(n5+i),1:nn)=X(1:(n5+i-1-t),1:nn);
      A56(t+1,1:nn)=A56(n5+i,1:nn);
      A56((t+2):(n5+i),1:nn)=Y(1:(n5+i-1-t),1:nn);
      U56(t+1,1:nn)=U56(n5+i,1:nn);
      U56((t+2):(n5+i),1:nn)=Z(1:(n5+i-1-t),1:nn);
      B56(t+1,1:n6)=B56(n5+i,1:n6);
      B56((t+2):(n5+i),1:n6)=W(1:(n5+i-1-t),1:n6);
      t=t+1;
   end
end

%All the orders of infinite zeros are 1.
if (n5==0)&(n6>0)
   U56(1:nn,1:nn)=inv(E56(1:nn,1:nn));
   V56(1:nn,1:nn)=eye(nn);
   A56(1:nn,1:nn)=inv(E56(1:nn,1:nn))*A56(1:nn,1:nn);
   CP56(1:nn,1:nn)=inv(C56(1:nn,1:nn));
   C56(1:nn,1:nn)=eye(nn);
   QB56(1:nn,1:nn)=inv(B56(1:nn,1:nn))*E56(1:nn,1:nn);
   B56(1:nn,1:nn)=eye(nn);
   E56(1:nn,1:nn)=eye(nn);
end

if size(p,2)>0
   for i=1:p(1)
      Ind(i)=Ind(i)+1;
   end
   for i=p(1)+1:n6
      Ind(i)=1;
   end
end
if n5==0
   Ind=ones(1,size(C56,1));
end