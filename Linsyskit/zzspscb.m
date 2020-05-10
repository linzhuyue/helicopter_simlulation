function [AA,BB,CC,Gs,Go,Gi,lv,rv,qv,dims,invGs]=zzspscb(A,B,C,tol,dc);

%ZZSPSCB  Decomposition of Strictly Proper System
%
%     [As,Bs,Cs,Gs,Go,Gi,lv,rv,qv,dims]=zzspscb(A,B,C[,tol,dc])
%
%     decomposes the system (A,B,C) into its scb structure (As,Bs,Cs).
%     
%     Inputs:
%        A, B, C  : state space matrices of a given system.       
%              dc : dc=0, for continuous-time system;
%                   dc=1, for discrete-time system.
%
%     Outputs:
%        As, Bs, Cs : state space matrices in an s.c.b.
%        Gs, Go, Gi : state, output and input transformations.
%        lv : = [l1, l2, ..., l_{pb}], observability index corresponding to x_b;
%        rv : = [r1, r2, ..., r_{mc}], controllability index corresponding to x_c;
%        qv : = [q1, q2, ..., q_{md}], infinite zeros, corresponding to x_d;
%        dims : = [nan, na0, nap, nb, nc, nd].
%
%     For details, see the manual and the references therein.
%
%     See also scb, scbraw.

if nargin==3
   tol=1e-8;dc=0;
elseif nargin==4
   dc=0;
end

%STEP SCB.1
[p,n]=size(C); [n,m]=size(B);
Ziall=reshape(C',1,n,p);
f=ones(p,1); Z=C; Zd=zeros(0,n); Wd=zeros(0,m); Zb=zeros(0,n);
%v=0; w=0; lv=0; rv=0; qv=0; alpha=ones(p,1);
v=0; w=0; lv=[]; rv=[]; qv=[]; alpha=ones(p,1);

%STEP SCB.2
while any(f)~=0
   for i=1:p
      if f(i)~=0
         Zi=Ziall(1:alpha(i),:,i);
         CiaiB=Zi(alpha(i),:)*B;
         if rank([Wd;CiaiB],tol)>rank(Wd,tol) %Case 1
            f(i)=0;
            Zd=[Zd;Zi]; Wd=[Wd;CiaiB];
            v=v+1; qv(v)=alpha(i);
            Zdj(1:size(Zi,1),:,v)=Zi;
         else %Case 2
            if isempty(Wd)
               vCi=Zi;
            else
               betai=CiaiB*pinv(Wd);
               for k=1:alpha(i)
                  t2=zeros(1,n);
                  for j=1:v
                     t1=qv(j)-alpha(i)+k;
                     if t1>0
                        t2=t2+betai(j)*Zdj(t1,:,j);
                     end
                  end
                  vCi(k,:)=Zi(k,:)-t2;
               end
            end
            vCiai1=vCi(alpha(i),:)*A;
            if rank([Z;vCiai1],tol)>rank(Z,tol) %Sub-case 2.1
               Z=[Z;vCiai1]; Zi=[vCi;vCiai1]; Ziall(1:1+alpha(i),:,i)=Zi;
               alpha(i)=alpha(i)+1;
            else %Sub-case 2.2
               f(i)=0; Zb=[Zb;vCi]; w=w+1; lv(w)=alpha(i);
            end
         end
      end
   end
end

%STEP SCB.3
md=v; nd=sum(qv); pb=w; nb=sum(lv); no=n-nd-nb;
tt=[Zb;Zd]; [t1,t1,Zo]=svd(tt); Zo=Zo(:,size(tt,1)+1:n);%Zo=null([Zb;Zd]);
S=[Zo';Zb;Zd];
[t1,t1,Wc]=svd(Wd); Wc=Wc(:,size(Wd,1)+1:m);%Wc=null(Wd);
W=[Wd;Wc'];
CCt=zeros(size(C)); BBt=zeros(size(B));
for i=1:md
   CCt(i,n-nd+sum(qv(1:i-1))+1)=1;
   BBt(n-nd+sum(qv(1:i)),i)=1;
end
for i=1:pb
   CCt(md+i,no+sum(lv(1:i-1))+1)=1;
end
Gs=inv(S);
invGs=S;
Go=C*inv(S)*CCt';
Gi=inv(W);

%STEP SCB.4
AA=invGs*A*Gs;
BB=invGs*B*Gi;
CC=inv(Go)*C*Gs;
%AA=inv(Gs)*A*Gs;
%BB=inv(Gs)*B*Gi;
%CC=inv(Go)*C*Gs;

%STEP SCB.5
t1=eye(n); t2=eye(p);
for i=1:pb
   r1=sum(lv(1:i-1))+no; r2=sum(lv(1:i))+no;
   for k=1:lv(i)
      for s=1:pb
         e1=sum(lv(1:s-1))+no;
         for j=lv(i)+2-k:min(lv(i)+1,lv(s))
            t1(r1+k,e1-lv(i)+j-1+k)=-AA(r2,e1+j);
         end
         if lv(s)>lv(i)
            t2(md+i,md+s)=-AA(r2,e1+lv(i)+1);
         end
      end      
      for s=1:md
         e1=sum(qv(1:s-1))+no+nb;
         for j=lv(i)+2-k:min(lv(i)+1,qv(s))
            t1(r1+k,e1-lv(i)+j-1+k)=-AA(r2,e1+j);
         end
         if qv(s)>lv(i)
            t2(md+i,s)=-AA(r2,e1+lv(i)+1);
         end
      end
   end
end
Gs=Gs*inv(t1);
invGs=t1*invGs;
Go=Go*inv(t2);
AA=invGs*A*Gs;
BB=invGs*B*Gi;
CC=inv(Go)*C*Gs;

%STEP SCB.6 and STEP.7
Bod=BB(1:no,1:md);Bd=BB(n-nd+1:n,1:md);
t1=eye(n);t1(1:no,n-nd+1:n)=Bod*Bd';
Gs=Gs*t1;
invGs=inv(t1)*invGs;
AA=invGs*A*Gs;
BB=invGs*B*Gi;
CC=inv(Go)*C*Gs;

%STEP SCB.8
t1=eye(n);t3=AA;
for i=0:max([qv lv])-2
   t2=eye(n);
   for s=1:pb
      if lv(s)>i+1
         t2(1:no,no+sum(lv(1:s))-i-1)=t3(1:no,no+sum(lv(1:s))-i);
      end
   end
   for s=1:md
      if qv(s)>i+1
         t2(1:no,no+nb+sum(qv(1:s))-i-1)=t3(1:no,no+nb+sum(qv(1:s))-i);
      end
   end
   t1=t1*t2;
   t3=inv(t1)*AA*t1;
end
Gs=Gs*t1;
invGs=inv(t1)*invGs;
AA=invGs*A*Gs;
BB=invGs*B*Gi;
CC=inv(Go)*C*Gs;

%STEP SCB.9
Aoo=AA(1:no,1:no);Boc=BB(1:no,md+1:m);
[t1,t2,Tos,Toi,na,rv]=csd(Aoo,Boc,tol);
nc=no-na;mc=m-md;
t1=eye(n);t1(1:no,1:no)=Tos;
t2=eye(m);t2(md+1:m,md+1:m)=Toi;
t3=eye(n);tt=t3(na+1:na+nc,:);
t3(na+1:na+nb,:)=t3(na+nc+1:na+nc+nb,:);
t3(na+nb+1:na+nb+nc,:)=tt;
Gs=Gs*t1*inv(t3);
invGs=t3*inv(t1)*invGs;
Gi=Gi*t2;
AA=invGs*A*Gs;
BB=invGs*B*Gi;
CC=inv(Go)*C*Gs;

%STEP SCB.10
nan=0; na0=0; nap=0;
if na>0
   Aa=AA(1:na,1:na);
   if dc==0
      [Aa,tt,nan,na0,nap]=ssd(Aa,tol);
   else
      [Aa,tt,nan,na0,nap]=dssd(Aa,tol);
   end
   t1=eye(n);t1(1:na,1:na)=tt;
   Gs=Gs*t1;
   invGs=inv(t1)*invGs;
   AA=invGs*A*Gs;
   BB=invGs*B*Gi;
   CC=inv(Go)*C*Gs;
end
dims=[nan,na0,nap,nb,nc,nd];

%Make the elements should be zeros be zeros.
[At,Bt,Ct,At1]=zzstandabc(na,lv,rv,qv);

AA=AA.*At;
AA=AA-AA.*At1+At1;
BB=BB.*Bt;CC=CC.*Ct;
err_SPSCB_A=norm(AA-invGs*A*Gs);
err_SPSCB_B=norm(BB-invGs*B*Gi);
err_SPSCB_C=norm(CC-inv(Go)*C*Gs);
%if (err_SPSCB_A>1e-8)|(err_SPSCB_B>1e-8)|(err_SPSCB_C>1e-8)
%   err_SPSCB_A,err_SPSCB_B,err_SPSCB_C
%end
%err_spscb=[err_SPSCB_A,err_SPSCB_B,err_SPSCB_C]