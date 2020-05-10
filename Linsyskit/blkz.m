function z=blkz(A,B,C,D,tol)

%BLKZ  Blocking Zeros of Multivariable Systems
%
%      bzero = blkz(A,B,C,D)
%
%      returns blocking zeros of a system characterized by (A,B,C,D).
%      A blocking zero, say 'alpha', is such that
%
%          H(alpha) = C (alpha*I - A)^{-1} B + D = 0
%
%      A blocking zero is also an invariant zero.
%
%      See also INVZ.

if isempty(A)|isempty(B)|isempty(C)|isempty(D)
   disp(' Sorry, there is empty matrix in A,B,C, ......')
   z=[];
   return
end   

if nargin==4,
   tol=1e-8;
end;
[ns,nu]=size(B);
tol1=10*ns*norm(A,1)*eps;
[Am,Bm,Cm,t,kk]=ctrbf(A,B,C,tol1);
kk=sum(kk);
Am=Am(ns-kk+1:ns,ns-kk+1:ns);
Bm=Bm(ns-kk+1:ns,:);
Cm=Cm(:,ns-kk+1:ns);
ns=kk;
if ns
   [Am,Bm,Cm,t,kk]=obsvf(Am,Bm,Cm,tol1);
   kk=sum(kk);
   Am=Am(ns-kk+1:ns,ns-kk+1:ns);
   Bm=Bm(ns-kk+1:ns,:);
   Cm=Cm(:,ns-kk+1:ns);
end
Dm=D;

dc=0;
[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=scb(Am,Bm,Cm,Dm,tol,dc);
%type=2;dc=0;d11_eye=1;
%[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=zzscbchu(Am,Bm,Cm,Dm,tol,type,dc,d11_eye);

na=sum(dims(1:3));
mu=length(qv);
z=[];
if na>0
   Aa=AA(1:na,1:na);
   [AA,T,struc,Jeig,err_of_JCF]=jcf(Aa,sqrt(tol));
   for k=1:size(Jeig,1)
      am=sum(struc(k,:));
      gm=sum(struc(k,:)>0);
      if gm==mu
         t=min((struc==0)*100000+struc);
         z=[z;ones(t,1)*(Jeig(k,1)+i*Jeig(k,2))];
      end
   end
end
