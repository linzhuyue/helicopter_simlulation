function [AA,BB,CC,V,Go,Gi,Ind,dims,U]=zzescbedz(E,A,B,C,type,tol)

% Structure Decomposition of Strictly Proper System (E,A,B,C,0)
%
% [AA,BB,CC,V,Go,Gi,Ind,dims,U]=zzescbedz(E,A,B,C,type,tol)

if nargin==4
   type=0; tol=1e-8;
elseif nargin==5,
   tol=1e-8;
end

n=size(A,1);
m=size(B,2); p=size(C,1);
if n==0
   AA=[];BB=[];CC=[];V=[];Go=[];Gi=[];Ind=[];dims=zeros(1,6);U=[];
   return;
end

%A-STEP.1
[U1,B6,mu6]=zzrowdown(B,tol);
[V1,C6,tao6]=zzrowdown(C',tol);V1=V1';C6=C6';
U=U1;V=V1;
if isempty(U)
   U=eye(n);
end
if isempty(V)
   V=eye(n);
end
EE=U*E*V;
AA=U*A*V;
E11=EE(1:n-mu6,1:n-tao6);
A11=AA(1:n-mu6,1:n-tao6);

%A-STEP.2
[E11,A11,U2,V2,nn,ll,na0]=zzfiguptri(E11,A11,tol,type);
U2I=eye(n);U2I(1:size(U2,1),1:size(U2,1))=U2;
V2I=eye(n);V2I(1:size(V2,1),1:size(V2,1))=V2;
U=U2I*U;V=V*V2I;
EE=U2I*EE*V2I;
AA=U2I*AA*V2I;

n1=ll(1);n2=ll(2);n3=ll(3);tao42=ll(4);n5=ll(5);
mu12=nn(1);n4=nn(4);

%A-STEP.3
E11=EE(1:mu12,1:n1);
A11=AA(1:mu12,1:n1);
E61=EE(sum(nn)+1:sum(nn)+mu6,1:n1);
A61=AA(sum(nn)+1:sum(nn)+mu6,1:n1);
U3=zzrowup([E11;E61]);
n6=mu12+mu6-n1;

U3I(1:n1,1:mu12)=U3(1:n1,1:mu12);
U3I(1:n1,n-mu6+1:n)=U3(1:n1,mu12+1:mu12+mu6);
U3I(n-n6+1:n,1:mu12)=U3(n1+1:n1+n6,1:mu12);
U3I(n-n6+1:n,n-mu6+1:n)=U3(n1+1:n1+n6,mu12+1:mu12+mu6);
U3I(n1+1:n-n6,mu12+1:n-mu6)=eye(n-n1-n6);
U=U3I*U;
EE=U3I*EE;
AA=U3I*AA;
BB=U*B;
B1=BB(1:n1,:);
B6=BB(n-n6+1:n,:);

E44=EE(n1+n2+n3+1:n1+n2+n3+n4,n1+n2+n3+1:n1+n2+n3+tao42);
E46=EE(n1+n2+n3+1:n1+n2+n3+n4,n1+n2+n3+tao42+n5+1:n1+n2+n3+tao42+n5+tao6);
V3=zzrowup([E44,E46]');V3=V3';
nt6=tao42+tao6-n4;
V3I=eye(n1+n2+n3);
V3I(n1+n2+n3+1:n1+n2+n3+tao42,n1+n2+n3+1:n1+n2+n3+n4)=V3(1:tao42,1:n4);
V3I(n1+n2+n3+1:n1+n2+n3+tao42,n-nt6+1:n)=V3(1:tao42,n4+1:n4+nt6);
V3I(n1+n2+n3+tao42+1:n-tao6,n1+n2+n3+n4+1:n-nt6)=eye(n5);
V3I(n-tao6+1:n,n1+n2+n3+1:n1+n2+n3+n4)=V3(tao42+1:tao42+tao6,1:n4);
V3I(n-tao6+1:n,n-nt6+1:n)=V3(tao42+1:tao42+tao6,n4+1:n4+nt6);
V=V*V3I;
EE=EE*V3I;
AA=AA*V3I;
CC=C*V;
C4=CC(:,n1+n2+n3+1:n1+n2+n3+n4);
C6=CC(:,n-nt6+1:n);

%A-STEP.4
W=zzrowup(B6');W=W';
P=zzrowup(C6);
BB=BB*W;
CC=P*CC;

%A-STEP.5
B11=BB(1:n1,1:n6);
B61=BB(n-n6+1:n,1:n6);
C14=CC(1:n6,n1+n2+n3+1:n1+n2+n3+n4);
C16=CC(1:n6,n-n6+1:n);
U5=zzrowdown([B11;B61]);
V5=zzrowdown([C14,C16]');V5=V5';
U5I(1:n1,1:n1)=U5(1:n1,1:n1); U5I(1:n1,n-n6+1:n)=U5(1:n1,n1+1:n1+n6);
U5I(n1+1:n-n6,n1+1:n-n6)=eye(n2+n3+n4+n5);
U5I(n-n6+1:n,1:n1)=U5(n1+1:n1+n6,1:n1); U5I(n-n6+1:n,n-n6+1:n)=U5(n1+1:n1+n6,n1+1:n1+n6);
V5I=eye(n1+n2+n3);
V5I(n1+n2+n3+1:n1+n2+n3+n4,n1+n2+n3+1:n1+n2+n3+n4)=V5(1:n4,1:n4);
V5I(n1+n2+n3+1:n1+n2+n3+n4,n-n6+1:n)=V5(1:n4,n4+1:n4+n6);
V5I(n-n5-n6+1:n-n6,n-n5-n6+1:n-n6)=eye(n5);
V5I(n-n6+1:n,n1+n2+n3+1:n1+n2+n3+n4)=V5(n4+1:n4+n6,1:n4);
V5I(n-n6+1:n,n-n6+1:n)=V5(n4+1:n4+n6,n4+1:n4+n6);
U=U5I*U;
V=V*V5I;
EE=U5I*EE*V5I;
AA=U5I*AA*V5I;
BB=U5I*BB;
CC=CC*V5I;
nn=[n1 n2 n3 n4 n5 n6];

%B-STEP.1
Et123=EE(1:n1+n2+n3,1:n1+n2+n3);
Et6=EE(n-n6+1:n,1:n1+n2+n3);
Utt=zzrowup([Et123;Et6]);
Ut=eye(n);
Ut(n-n6+1:n,1:n1+n2+n3)=Utt(n1+n2+n3+1:n1+n2+n3+n6,1:n1+n2+n3);
Ut(n-n6+1:n,n-n6+1:n)=Utt(n1+n2+n3+1:n1+n2+n3+n6,n1+n2+n3+1:n1+n2+n3+n6);
E44=EE(n1+n2+n3+1:n1+n2+n3+n4,n1+n2+n3+1:n1+n2+n3+n4);
E46=EE(n1+n2+n3+1:n1+n2+n3+n4,n-n6+1:n);
Vtt=zzrowup([E44,E46]');Vtt=Vtt';
Vt=eye(n);
Vt(n1+n2+n3+1:n1+n2+n3+n4,n-n6+1:n)=Vtt(1:n4,n4+1:n4+n6);
Vt(n-n6+1:n,n-n6+1:n)=Vtt(n4+1:n4+n6,n4+1:n4+n6);
U=Ut*U;V=V*Vt;
EE=Ut*EE*Vt;
AA=Ut*AA*Vt;
BB=Ut*BB;
CC=CC*Vt;

%B-STEP.2
B6162t=BB(n-n6+1:n,:);
C1626t=CC(:,n-n6+1:n);
Wt=zzrowup(B6162t');Wt=Wt';
Wt(:,1:n6)=[eye(n6);zeros(m-n6,n6)];
Pt=zzrowup(C1626t);
Pt(1:n6,:)=[eye(n6),zeros(n6,p-n6)];
Zb=W*Wt;Zc=Pt*P;
BB=BB*Wt;
CC=Pt*CC;
B12=BB(1:n1,n6+1:m);
C24=CC(n6+1:p,n1+n2+n3+1:n1+n2+n3+n4);

%B-STEP.3-----refine by using Generalized Sylvester equations
E11=EE(1:n1,1:n1);
E22=EE(n1+1:n1+n2,n1+1:n1+n2);
E33=EE(n1+n2+1:n1+n2+n3,n1+n2+1:n1+n2+n3);
E44=EE(n1+n2+n3+1:n1+n2+n3+n4,n1+n2+n3+1:n1+n2+n3+n4);
A11=AA(1:n1,1:n1);
A22=AA(n1+1:n1+n2,n1+1:n1+n2);
A33=AA(n1+n2+1:n1+n2+n3,n1+n2+1:n1+n2+n3);
A44=AA(n1+n2+n3+1:n1+n2+n3+n4,n1+n2+n3+1:n1+n2+n3+n4);
t1=eig(inv(E11)*A11);
t2=eig(inv(E22)*A22);
t3=eig(inv(E33)*A33);
eit1=[t1;t2;t3];
eit2=[t2;t3];


X1I=eye(n);X2I=eye(n);X3I=eye(n);X4I=eye(n);
Y1I=eye(n);Y2I=eye(n);Y3I=eye(n);Y4I=eye(n);
   
if (n1+n2+n3~=0)&(n4+n5+n6~=0)
   Et=EE(1:n1+n2+n3,1:n1+n2+n3);
   et1=EE(n1+n2+n3+1:n1+n2+n3+n4+n5,n1+n2+n3+1:n);
   et1=[et1;zeros(n6,n4+n5+n6)];
   et2=-EE(1:n1+n2+n3,n1+n2+n3+1:n);
   At=AA(1:n1+n2+n3,1:n1+n2+n3);
   at2=-AA(1:n1+n2+n3,n1+n2+n3+1:n);
   
   if (norm(A44)>tol)&(norm(C24)>tol), t1=norm(A44)/norm(C24); else t1=1; end
   t1=t1*5;
   eigt1=zeros(n1+n2+n3,0);tm=0;
   for k=1:n4, eigt1=[eigt1,eit1]; end
   errm=Inf;
   for kk=1:20
      tt=0;
      while tt<tol*1000
         N=2*t1*rand(n4,p-n6)-t1;
         t3=eig(inv(E44*(A44-N*C24))); eigt2=zeros(0,n4);
         for k=1:n1+n2+n3, eigt2=[eigt2;t3']; end
         tt=min(min(abs(eigt2-eigt1)));
      end
      at1=AA(n1+n2+n3+1:n1+n2+n3+n4+n5,n1+n2+n3+1:n);
      at1(1:n4,1:n4)=at1(1:n4,1:n4)-N*C24;
      at1=[at1;zeros(n6,n4+n5),-C16];
      [Y1t,x1kt,err]=zzgensyl(Et,et1,et2,At,at1,at2);
      if err<errm
         errm=err;
         Y1=Y1t;
         x1k=x1kt;
      end
   end   
   X1=x1k(:,1:n4+n5);
   X1I(1:n1+n2+n3,n1+n2+n3+1:n-n6)=X1;
   Y1I(1:n1+n2+n3,n1+n2+n3+1:n)=Y1;
   AA=X1I*AA*Y1I;
   EE=X1I*EE*Y1I;
   EE(1:n1+n2+n3,n1+n2+n3+1:n)=0;
   AA(1:n1+n2+n3,n1+n2+n3+n4+1:n-n6)=0;
end

if (n4~=0)&(n5+n6~=0)
   Et=EE(n1+n2+n3+1:n1+n2+n3+n4,n1+n2+n3+1:n1+n2+n3+n4);
   et1=EE(n1+n2+n3+n4+1:n,n1+n2+n3+n4+1:n-n6);
   et1=[et1,zeros(n5+n6,size(B61,2))];
   et2=-EE(n1+n2+n3+n4+1:n,n1+n2+n3+1:n1+n2+n3+n4);
   At=AA(n1+n2+n3+1:n1+n2+n3+n4,n1+n2+n3+1:n1+n2+n3+n4);
   at1=AA(n1+n2+n3+n4+1:n,n1+n2+n3+n4+1:n-n6);
   at1=[at1,[zeros(n5,size(B61,2));-B61]];
   at2=-AA(n1+n2+n3+n4+1:n,n1+n2+n3+1:n1+n2+n3+n4);
   [Y2f,X2]=zzgensyl(et1,Et,et2,at1,At,at2);
   Y2=Y2f(1:n5,:);
   Y2I(n1+n2+n3+n4+1:n-n6,n1+n2+n3+1:n1+n2+n3+n4)=Y2;
   X2I(n1+n2+n3+n4+1:n,n1+n2+n3+1:n1+n2+n3+n4)=X2;
   AA=X2I*AA*Y2I;
   EE=X2I*EE*Y2I;
   EE(n1+n2+n3+n4+1:n,n1+n2+n3+1:n1+n2+n3+n4)=0;
   AA(n1+n2+n3+n4+1:n-n6,n1+n2+n3+1:n1+n2+n3+n4)=0;
end;

if (n1~=0)&(n2+n3~=0)
   Et=EE(n1+1:n1+n2+n3,n1+1:n1+n2+n3);
   et1=EE(1:n1,1:n1);
   et2=-EE(1:n1,n1+1:n1+n2+n3);
   At=AA(n1+1:n1+n2+n3,n1+1:n1+n2+n3);
   at2=-AA(1:n1,n1+1:n1+n2+n3);
   
   if (norm(A11)>tol)&(norm(B12)>tol), t1=norm(A11)/norm(B12); else t1=1; end
   t1=t1*5;
   eigt1=zeros(n2+n3,0);tm=0;
   for k=1:n1, eigt1=[eigt1,eit2]; end
   
   errm=Inf;
   for i=1:20
      tt=0;
      while tt<tol*1000
         M=2*t1*rand(m-n6,n1)-t1;
         t3=eig(inv(E11)*(A11-B12*M)); eigt2=zeros(0,n1);
         for k=1:n2+n3, eigt2=[eigt2;t3']; end
         tt=min(min(abs(eigt2-eigt1)));
      end
      at1=AA(1:n1,1:n1)-B12*M;
      [Y3t,X3t,err]=zzgensyl(et1,Et,et2,at1,At,at2);      
      if err<errm
         errm=err;
         Y3=Y3t;
         X3=X3t;
      end
   end   
   Y3I=eye(n);
   Y3I(1:n1,n1+1:n1+n2+n3)=Y3;
   X3I=eye(n);
   X3I(1:n1,n1+1:n1+n2+n3)=X3;
   AA=X3I*AA*Y3I;
   EE=X3I*EE*Y3I;
   EE(1:n1,n1+1:n1+n2+n3)=0;
end


if 0&(n2~=0)&(n3~=0)
   Et=EE(n1+1:n1+n2,n1+1:n1+n2);
   et1=EE(n1+n2+1:n1+n2+n3,n1+n2+1:n1+n2+n3);
   et2=-EE(n1+1:n1+n2,n1+n2+1:n1+n2+n3);
   At=AA(n1+1:n1+n2,n1+1:n1+n2);
   at1=AA(n1+n2+1:n1+n2+n3,n1+n2+1:n1+n2+n3);
   at2=-AA(n1+1:n1+n2,n1+n2+1:n1+n2+n3);
   [Y4,X4]=zzgensyl(Et,et1,et2,At,at1,at2);
   Y4I=eye(n);
   Y4I(n1+1:n1+n2,n1+n2+1:n1+n2+n3)=Y4;
   X4I=eye(n);
   X4I(n1+1:n1+n2,n1+n2+1:n1+n2+n3)=X4;
   AA=X4I*AA*Y4I;
   EE=X4I*EE*Y4I;
   AA(n1+1:n1+n2,n1+n2+1:n1+n2+n3)=0;
   EE(n1+1:n1+n2,n1+n2+1:n1+n2+n3)=0;
end;

U=X4I*X3I*X2I*X1I*U;
V=V*Y1I*Y2I*Y3I*Y4I;   

%Transform the state of system in the order of (xa xb xc xd), and make the elements should be zeros be zeros.
dims=[n2,na0,n3-na0,n4,n1,n5+n6];
t1=eye(n);
UU=[t1(n1+1:n1+n2+n3+n4,:);t1(1:n1,:);t1(n-n5-n6+1:n,:)];
U=UU*U;V=V*UU';
EE=U*E*V;
AA=U*A*V;
BB=U*B*Zb;
CC=Zc*C*V;
na=n2+n3;nb=n4;nc=n1;nd=n5+n6;
Et=blkdiag(ones(na),ones(nb),ones(nc),ones(nd));
At=ones(n);
At(na+1:na+nb,1:na)=0;
At(1:na+nb,na+nb+1:na+nb+nc)=0;
Bt=zeros(n,m);
Bt(n-nd+1:n,1:n6)=1;
Bt(na+nb+1:n-nd,n6+1:m)=1;
Ct=zeros(p,n);
Ct(1:n6,n-nd+1:n)=1;
Ct(n6+1:p,na+1:na+nb)=1;

%%%
At(n-n5-n6+1:n-n6,1:n1+n2+n3+n4)=0;
At(1:n1+n2+n3+n4,n-n5-n6+1:n-n6)=0;
At(n2+1:n2+n3,1:n2)=0;
%At(1:n2,n2+1:n2+n3)=0;%!!!!!!!
Bt(n-n6-n5+1:n-n6,:)=0;
Ct(:,n-n6-n5+1:n-n6)=0;
%%%

EE=Et.*EE;
AA=At.*AA;
if ~isempty(BB)
   BB=Bt.*BB;
end
if ~isempty(CC)
   CC=Ct.*CC;
end

%Extraction of Infinite Zero Structure
Ind=[];
if (type>=1) & (n5+n6>0)
   E56=EE(n-n5-n6+1:n,n-n5-n6+1:n);
   A56=AA(n-n5-n6+1:n,n-n5-n6+1:n);
   B56=BB(n-n5-n6+1:n,1:n6);
   C56=CC(1:n6,n-n5-n6+1:n);
   [U56,V56,QB56,CP56,A56,B56,C56,Ind]=zzinfzerogsvd(E56,A56,B56,C56,n5,n6,tol);
   U5I=eye(n);V5I=eye(n);
   U5I(n-n5-n6+1:n,n-n5-n6+1:n)=U56;
   V5I(n-n5-n6+1:n,n-n5-n6+1:n)=V56;
   EE=U5I*EE*V5I;AA=U5I*AA*V5I;
   EE(n-n5-n6+1:n,n-n5-n6+1:n)=eye(n5+n6);;
   AA(n-n5-n6+1:n,n-n5-n6+1:n)=A56;
   BB(n-n5-n6+1:n,1:n6)=B56;
   CC(1:n6,n-n5-n6+1:n)=C56;
   U=U5I*U;V=V*V5I;
   Zb(1:m,1:n6)=Zb(1:m,1:n6)*QB56;
   Zc(1:n6,1:p)=CP56*Zc(1:n6,1:p);
end

V=V*inv(EE);
AA=AA*inv(EE);
CC=CC*inv(EE);
EE=eye(n);

invGo=Zc;
Go=inv(Zc);
Gi=Zb;

err_escbedz_A=norm(U*A*V-AA);err_escbedz_B=norm(U*B*Gi-BB);err_escbedz_C=norm(inv(Go)*C*V-CC);
%if any([err_escbedz_A,err_escbedz_B,err_escbedz_C]>1e-8)
%   err_escbedz_A,err_escbedz_B,err_escbedz_C
%end
%err_escbedz=[err_escbedz_A,err_escbedz_B,err_escbedz_C]