function [A,B,C,At]=zzstandabc(na,bv,cv,dv)

nb=sum(bv);nc=sum(cv);nd=sum(dv);
pb=sum(bv~=0);mc=sum(cv~=0);md=sum(dv~=0);

Aa=ones(na,na);

Abt=[];Act=[];Adt=[];

Ab=diag(ones(1,nb-1),1);Cb=zeros(pb,nb);
for i=1:pb
   Cb(i,sum(bv(1:i-1))+1)=1;
   if i~=pb Ab(sum(bv(1:i)),sum(bv(1:i))+1)=0; end
end
if pb~=0
   Abt=Ab;
end
Ab=Ab+ones(nb,pb)*Cb;

Ac=diag(ones(1,nc-1),1);Bc=zeros(nc,mc);
for i=1:mc
   Bc(sum(cv(1:i)),i)=1;
   if i~=mc Ac(sum(cv(1:i)),sum(cv(1:i))+1)=0; end
end
if mc~=0
   Act=Ac;
end
Ac=Ac+Bc*ones(mc,nc);

Ad=diag(ones(1,nd-1),1);Cd=zeros(md,nd);Bd=zeros(nd,md);
for i=1:md
   Cd(i,sum(dv(1:i-1))+1)=1;Bd(sum(dv(1:i)),i)=1;
   if i~=md Ad(sum(dv(1:i)),sum(dv(1:i))+1)=0; end
end
if md~=0
   Adt=Ad;
end
Ad=Ad+ones(nd,md)*Cd+Bd*ones(md,nd)-Bd*Cd;

LabCb=ones(na,pb)*Cb;LcbCb=ones(nc,pb)*Cb;LCd=ones(na+nb+nc,md)*Cd;
BcEca=Bc*ones(mc,na);BdEd=Bd*ones(md,na+nb+nc);

At=blkdiag(zeros(na,na),Abt,Act,Adt);
A=[[Aa;zeros(nb,na);BcEca],[LabCb;Ab;LcbCb],[zeros(na+nb,nc);Ac],LCd;[BdEd,Ad]];
A=(A~=0);
B=[zeros(na+nb,mc+md);zeros(nc,md),Bc;Bd,zeros(nd,mc)];
C=[zeros(md,na+nb+nc),Cd;zeros(pb,na),Cb,zeros(pb,nc+nd)];
D=zeros(pb+md,mc+md);