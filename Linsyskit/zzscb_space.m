function S = zzscb_space(A,B,C,D,flag,dc,tol,lambda)

%SCB_SPACE  Geometric Subspace from SCB
%
%           S = zzscb_space(A,B,C,D,flag,dc,tol,lambda)
%
%           Including: v_star v_minus v_plus s_star s_minus s_plus r_star n_star s_lambda and v_lambda.

if nargin==5,
   dc=0;tol=1e-8;lambda=0;
elseif nargin==6
   tol=1e-8;lambda=0;
elseif nargin==7
   lambda=0;
end;

%[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=scb(A,B,C,D,tol,dc);

type=0;d11_eye=0;
[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);

n=size(A,1);
nan=dims(1); na0=dims(2); nap=dims(3);
na=sum(dims(1:3)); nb=dims(4); nc=dims(5); nd=dims(6);

if na>0 & dc==0 
   [D,T,nan,na0,nap] = ssd(AA(1:na,1:na),tol);
   T(na+1:n,na+1:n)=eye(n-na);
   Gs=Gs*T;
   AA=inv(T)*AA*T;
   BB=inv(T)*BB;
   CC=CC*T;
end

if na>0 & dc==1 
   [D,T,nan,na0,nap] = dssd(AA(1:na,1:na),tol);
   T(na+1:n,na+1:n)=eye(n-na);
   Gs=Gs*T;
   AA=inv(T)*AA*T;
   BB=inv(T)*BB;
   CC=CC*T;
end

if flag=='v_star'
   t=zeros(n,na+nc); t(1:na,1:na)=eye(na); t(na+nb+1:na+nb+nc,na+1:na+nc)=eye(nc);
   S=Gs*t;
   return;   
elseif flag=='v_minu'
   t=zeros(n,nan+na0+nc); t(1:nan+na0,1:nan+na0)=eye(nan+na0); t(na+nb+1:na+nb+nc,nan+na0+1:nan+na0+nc)=eye(nc);
   S=Gs*t;
   return;
elseif flag=='v_plus'
   t=zeros(n,nap+nc); t(nan+na0+1:na,1:nap)=eye(nap); t(na+nb+1:na+nb+nc,nap+1:nap+nc)=eye(nc);
   S=Gs*t;
   return;
elseif flag=='s_star'
   t=zeros(n,nc+nd); t(na+nb+1:n,1:nc+nd)=eye(nc+nd);
   S=Gs*t;
   return;
elseif flag=='s_minu'
   t=zeros(n,nap+nc+nd); t(nan+na0+1:na,1:nap)=eye(nap); t(na+nb+1:n,nap+1:nap+nc+nd)=eye(nc+nd);
   S=Gs*t;
   return;
elseif flag=='s_plus'
   t=zeros(n,nan+na0+nc+nd); t(1:nan+na0,1:nan+na0)=eye(nan+na0); t(na+nb+1:n,nan+na0+1:nan+na0+nc+nd)=eye(nc+nd);
   S=Gs*t;
   return;
elseif flag=='r_star'
   t=zeros(n,nc); t(na+nb+1:n-nd,1:nc)=eye(nc);
   S=Gs*t;
   return;
elseif flag=='n_star'
   t=zeros(n,na+nc+nd); t(1:na,1:na)=eye(na); t(na+nb+1:n,na+1:na+nc+nd)=eye(nc+nd);
   S=Gs*t;
   return;
elseif flag=='s_lamb'
   md=length(qv);p=size(C,1);m=size(B,2);
   pb=p-md-m0;
   Yblambda=[];
   if nb~=0
      Ab=AA(na+1:na+nb,na+1:na+nb);Cb=CC(size(CC,1)-pb+1:size(CC,1),na+1:na+nb);
      tm=0;
      for k=1:20
         Kb=rand(nb,pb);
         tt=min(abs(eig(Ab-Kb*Cb)-lambda));
         if tt>tm
            Kbm=Kb;tm=tt;
         end
      end
      Yblambda=null(Cb*inv(Ab+Kbm*Cb-lambda*eye(nb)));
   end
   t=blkdiag(lambda*eye(na)-AA(1:na,1:na),Yblambda,eye(nc+nd));
   S=Gs*t;
   return;
elseif flag=='v_lamb'
   md=length(qv);p=size(C,1);m=size(B,2);
   mc=m-md-m0;
   Xalambda=zznulltol(lambda*eye(na)-AA(1:na,1:na),tol);
   Xclambda=[];
   if nc~=0
      Acc=AA(na+nb+1:na+nb+nc,na+nb+1:na+nb+nc);Bc=BB(na+nb+1:na+nb+nc,m-mc+1:m);
      tm=0;
      for k=1:20
         Fc=rand(mc,nc);
         tt=min(abs(eig(Acc-Bc*Fc)-lambda));
         if tt>tm
            Fcm=Fc;tm=tt;
         end
      end
      Xclambda=inv(Acc+Bc*Fcm-lambda*eye(nc))*Bc;
   end
   t=zeros(n,size(Xalambda,2)+size(Xclambda,2));
   t(1:size(Xalambda,1),1:size(Xalambda,2))=Xalambda;
   t(size(Xalambda,1)+nb+1:size(Xalambda,1)+nb+size(Xclambda,1),size(Xalambda,2)+1:size(Xalambda,2)+size(Xclambda,2))=Xclambda;
   S=Gs*t;
   return;
end