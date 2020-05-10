function [AA,P,stru,flag]=zzrepvalueJ(A,tol)

%[J,P]=zzrepvalueJ(A[,tol])
%
%     transforms a square matrix A with one repeated real eigenvalue
%     into the Jordan form, i.e., J = inv(P)*A*P is in the Jordan form.
%

%      stru is a row vector containing the sizes of Jordan blocks.

%      Refer to the following for detail:
%
%      Computation of Generalized Eigenvectros, S. Bingulac & D.W. Luse
%      Computers & Elect. Engng, Vol. 15, No.1, pp.29-32, 1989

%     flag=1, the decomposition is reasonable for the given tol.
%     flag=0, the decomposition can not be carried out for the given tol.

if nargin==1,
   tol=1e-8;
end

flag=1;
n=size(A,1);
if isempty(A)
   AA=[];P=[];stru=[];
   return
end

eigen=mean(real(eig(A,'nobalance')));
B=A-eigen*eye(n);

P=[];
R=zzorthtol(B,tol);
qj=n-size(R,2);
ii=0;
N=eye(n);
stru=[];
qh=0;
while qh<qj
   ii=ii+1;
   N=B*N;
   Q=zznulltol(N,tol);
   pQ=size(Q,2);
   s=svd(R);
   pR=sum(s>tol);
   s=svd([R Q]);
   pRQ=sum(s>tol);
   ki=pRQ-pR;
   if (ki~=0)&(pQ~=0)
      uu=nchoosek(1:pQ,ki);
   else
      uu=[];
   end
   sw=inf;M=[];
   for kk=1:size(uu,1)
      Mt=[];
      for kkk=1:size(uu,2)
         Mt=[Mt,Q(:,uu(kk,kkk))];
      end
      s=cond([R Mt]);
      if s<sw
         M=Mt; sw=s;
      end
   end
   if ~isempty(M)
      for kk=1:size(M,2)
         Pl=[];
         mm=eye(n);
         for k=1:ii
            Pl=[mm*M(:,kk) Pl];
            mm=B*mm;
         end
         P=[P Pl];
         stru=[stru,size(Pl,2)];
      end
      R=[R M];
      qh=qh+size(M,2);
   end
end

[x,y]=size(P);
if (x~=n) | (y~=n) | (cond(P)>1e15)
   AA=A;P=eye(n);stru=[];flag=0;
   return;
end
AA=inv(P)*A*P;