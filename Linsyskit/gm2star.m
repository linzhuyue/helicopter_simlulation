function gms2 = gm2star(A,B,C,D,E,tol)

% gms2 = gm2star(A,B,C2,D2,E)
%
% Calculates the infimum or the best achievable performance of the H2 
% suboptimal control problem for the system:
%   .
%   x =  A x +  B u + E w
%   y =  x
%   h = C2 x + D2 u
%
% under all possible stabilizing state feedback controllers.
%
% See also gm2sos, h2state, h2out, gm8star and dgm2star.

% Modified by Ben M. Chen on April 17, 2020 at CUHK

if nargin==5
   tol=1e-8;
end

[A,B,C,D,Gs,Go,Gi,dims,lv,rv,qv,m0]=scb(A,B,C,D,tol);

%type=2;dc=0;d11_eye=1;
%[A,B,C,D,Gs,Go,Gi,dims,lv,rv,qv,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);

nan=dims(1);na0=dims(2);nap=dims(3);na=sum(dims(1:3));
nb=dims(4);nc=dims(5);nd=dims(6);
pb=length(lv);mc=length(rv);md=length(qv);
n=size(A,1);
p=size(C,1);
m=size(B,2);


gms2=[];
if na0~=0
   [n,m] = size(B);
   In = eye(n);
   Im = eye(m);
   Ip = eye(p);

   flag = 1;
   epsilon = 10;
   P = In;
   Q = In;
   gm2ss = 666;

   while flag == 1
       gm2s = gm2ss;
       epsilon = epsilon/10;
       C2t = [C; epsilon*In; zeros(m,n)];
       D2t = [D; zeros(n,m); epsilon*Im];
       
       P = h2care(A,B,C2t,D2t);
       evP = eig(P);
       if min(evP) < 0
           flag = 0;
       end

      gm2ss=sqrt(trace(E'*P*E));
  
      if abs(gm2s-gm2ss) < 0.01
          flag = 0;
      end
   end
   gms2 = gm2ss;
   return
end

if (nap+nb)==0
   gms2=0;
   return;
end

Ass=A(nan+1:n-nc-nd,nan+1:n-nc-nd);
B0s=B(nan+1:n-nc-nd,1:m0);
Cd=C(m0+1:m0+md,n-nd+1:n);
Lsd=A(nan+1:n-nc-nd,n-nd+1:n)*Cd';
Bs=[B0s Lsd];

Cs=C(:,nan+1:na+nb);Cs(1:m0,:)=0;
Cs=Go*Cs;
Ds=Go*[eye(m0+md);zeros(p-m0-md,m0+md)];

tE=inv(Gs)*E;
Es=tE(nan+1:nan+nap+nb,:);
[Ps,t,Fs,t] = care(Ass,Bs,Cs'*Cs,Ds'*Ds,Cs'*Ds,eye(size(Cs,2)));
et=real(eig(Ps));
if all(et>0)
   gms2=sqrt(trace(Es'*Ps*Es));
end
