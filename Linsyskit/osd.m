function [AA,CC,Ts,To,no,Oidx,err_OSD_A,err_OSD_C,invTs]=osd(A,C,tol)

%OSD  Observability Structural Decomposition for Unforced Systems
%
%     [At,Ct,Ts,To,uom,Oidx] = osd(A,C)
%
%     returns an observability structural decomposition for (A,C).
%
%     Input Parameters:
%            .
%            x = A x,   y = C x
%
%     Output Parameters:
%                           [ Ao *     0     ...  *    0     ]
%                           | 0  * I_{k_1-1} ...  *    0     |
%       At = inv(Ts)*A*Ts = | 0  *     0     ...  *    0     |
%                           | :  :     :      .   :    :     |
%                           | 0  *     0     ...  * I_{k_p-1}|
%                           [ 0  *     0     ...  *    0     ]
%
%                           [ 0  1     0     ...  0    0     ]
%       Ct = inv(To)*C*Ts = | :  :     :      .   :    :     |
%                           [ 0  0     0     ...  1    0     ]
%
%     uom = unobservable modes & Oidx = observability index of (C,A)
%
%     See OBVIDX, BDOSD and CSD.

if nargin==2
   tol=1e-8;
end

CCC=C;
pp=size(C,1);
n=size(A,1);

%check if C is full row rank.
[U,t,p]=zzrowup(C,tol);
if p==pp
   U=eye(pp);
else
   C=t;
end

if (p==0)|(n==0)
    AA=A;CC=CCC;Ts=eye(n);To=eye(pp);no=n;Oidx=[];
    err_OSD_A=0;err_OSD_C=0;invTs=eye(n);
    return;
end

%STEP OSD.1
[p,n]=size(C);
Ziall=reshape(C',1,n,p);
f=ones(p,1);
Z=C;
tZ=zeros(0,n);
w=0;
alpha=ones(p,1);
Oidx=[];
CCt=0*C;

%STEP OSD.2
while any(f)~=0
   for i=1:p
      if f(i)~=0
         Zi=Ziall(1:alpha(i),:,i);
         CiaiA=Zi(alpha(i),:)*A;
         if rank([Z;CiaiA],tol)>rank(Z,tol) %Case 1
            Z=[Z;CiaiA]; Zi=[Zi;CiaiA]; Ziall(1:1+alpha(i),:,i)=Zi;
            alpha(i)=alpha(i)+1;
         else %Case 2
            f(i)=0; tZ=[tZ;Zi]; w=w+1; Oidx(w)=alpha(i);
            To(w,i)=1; CCt(w,sum(Oidx(1:w-1))+1)=1;
         end
      end
   end
end

%STEP OSD.3
nob=size(tZ,1);no=n-nob;
[tt,tt,So]=svd(tZ);So=So(:,nob+1:n);%So=null(tZ)
S=[So';tZ];
CCt=[zeros(p,no),CCt(:,1:nob)];
Ts=inv(S);
invTs=S;
%AA=inv(Ts)*A*Ts;
AA=invTs*A*Ts;
%CC=inv(To)*C*Ts;

%STEP OSD.4
t1=eye(n);t2=eye(p);
for i=1:p
   r1=sum(Oidx(1:i-1))+no; r2=sum(Oidx(1:i))+no;
   for kk=1:Oidx(i)
      for s=1:p
         e1=sum(Oidx(1:s-1))+no;
         for j=Oidx(i)+2-kk:Oidx(s)
            t1(r1+kk,e1-Oidx(i)+j-1+kk)=-AA(r2,e1+j);
         end
      end
   end
end
Ts=Ts*inv(t1);
invTs=t1*invTs;
AA=invTs*A*Ts;
%CC=inv(To)*C*Ts;

%STEP OSD.5
t1=eye(n);t3=AA;
for i=1:p
   for j=Oidx(i)-1:-1:1
      t2=eye(n);
      t2(1:no,no+sum(Oidx(1:i-1))+j)=t3(1:no,no+sum(Oidx(1:i-1))+j+1);
      t1=t1*t2;
      t3=inv(t1)*AA*t1;
   end
end
Ts=Ts*t1;
invTs=inv(t1)*invTs;
%AA=inv(Ts)*A*Ts;
AA=invTs*A*Ts;
CC=inv(To)*C*Ts;

t1=CC*CCt';
To=To*t1;
%AA=inv(Ts)*A*Ts;
AA=invTs*A*Ts;
CC=inv(To)*C*Ts;

tt=zeros(size(AA));tt(1:no,1:no)=ones(no,no);
ttt=zeros(size(C));
if no~=n
   for i=1:length(Oidx)
      t1=no+sum(Oidx(1:i-1));t2=no+sum(Oidx(1:i));
      tt(:,t1+1)=ones(n,1);
      tt(t1+1:t2-1,t1+2:t2)=eye(Oidx(i)-1);
      ttt(i,t1+1)=1;
   end
end
AA=AA.*tt;
CC=CC.*ttt;

%C is not full row rank
CC=[CC;zeros(pp-p,n)];
To=U'*blkdiag(To,eye(pp-p));%To=U'*To;

err_OSD_A=norm(AA-invTs*A*Ts);
err_OSD_C=norm(CC-inv(To)*CCC*Ts);
%if err_OSD_A>1e-8 | err_OSD_C>1e-8
%   err_OSD_A,err_OSD_C
%end
%err_osd=[err_OSD_A,err_OSD_C]
