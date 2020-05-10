function [AA,BB,Ts,Ti,ks,err_of_BDCSD_A,err_of_BDCSD_B,flag] = bdcsd(A,B,tol,for_sd_ds)

%BDCSD  Block Diagonal Controllable Structural Decomposition
%
%       [At,Bt,Ts,Ti,ks] = bdcsd(A,B)
%
%       transforms a controllable pair (A,B) into the block diagonal
%       controllable structural decomposition form.
%
%       Input Parameters:
%               .
%               x = A x + B u
%
%       Output Parameters:
%
%                                [ A_1   0    0  ...   0  ]
%                                |  0   A_2   0  ...   0  |
%            At = inv(Ts)*A*Ts = |  0    0   A_3 ...   0  |
%                                |  :    :    :   .    :  |
%                                [  0    0    0  ...  A_k ]
%
%                                [ B_1   *    *  ...   *   * ]
%                                |  0   B_2   *  ...   *   * |
%            Bt = inv(Ts)*B*To = |  0    0   B_3 ...   *   * |
%                                |  :    :    :   .    :   : |
%                                [  0    0    0  ...  B_k  * ]
%       where
%                         [ 0  1  0 ... 0 ]         [ 0 ]
%                         | 0  0  1 ... 0 |         | 0 |
%                   A_i = | :  :  :  .  : |,  B_i = | : |
%                         | 0  0  0 ... 1 |         | 0 |
%                         [ *  *  * ... * ]         [ 1 ]
%
%            ks: contains the sizes of the blocks A_1, A_2, ..., A_k.
%
%       See also CSD and BDOSD.

%       flag=0, the system is not completely controllable,
%       flag=1, the system is completely controllable.

if nargin==2
   tol=1e-8;
end

BBB=B;
n=size(A,1);
flag=1;

[t,t,t,t,nox]=ctrbf(A,B,zeros(1,n),tol);
nox=n-sum(nox);
if nox>0
   disp(' ')
   disp('  ***********************************************************************************')
   disp('  *** The system is not completely controllable or is not completely observable ***') 
   disp('  ***********************************************************************************') 
   AA=A;BB=B;
   Ts=eye(size(A,1));Ti=eye(size(B,2));ks=[];
   err_of_BDCSD_A=0;err_of_BDCSD_B=0;flag=0;
   return
end

mm=size(B,2);

%check if B is full column rank.
[U,t,m]=zzrowup(B',tol);t=t';U=U';

if m==mm
   U=eye(mm);
else
   B=t;
%   disp('  *********************************')
%   disp('  *** B is not full column rank ***')
%   disp('  *********************************')
end

if (m==0)|(n==0)
   AA=A;BB=B;Ts=eye(size(A,1));Ti=eye(mm);ks=0;
   err_of_BDCSD_A=0;err_of_BDCSD_B=0; return;
end

a=A;
b=B;
fla=1;

ks=[];
Ts=eye(n);
invTs=eye(n);
Ti=eye(m);
while fla==1
   if nargin==4
      [At,Q,xindex]=rjd(a,1234567,1);
   else
      [At,Q,xindex]=rjd(a,tol,1);
   end
   tB=inv(Q)*b;
   [ell,k]=size(xindex);
   xsize=sum(xindex,2);
   contr=sum(xindex(:,1)); 
   [ni,mi]=size(tB);
   t1=inf; k1=0;
   while (k1<10) | (t1==inf)
      k1=k1+1;
      T2=rand(mi,1); T2(1)=1;
      t2=svd(ctrb(a,b*T2));
      if t2(contr)>1e-7
         t2=t2(1)/t2(contr);
         if t1>t2
            t1=t2; T1=T2;
         end
      end
   end
   tB1=tB*T1;
   Toi=[T1 [zeros(1,mi-1); eye(mi-1)]]; 
   Tssi=[];
   invTssi=[];
   TT=eye(ni);
   P1=[];P2=[]; 
   for s=1:ell 
      dimx=sum(xsize(1:s));
      dimx_1=sum(xsize(1:s-1));
      dims=xsize(s); 
      slambdai1=xindex(s,1);
      at=At(dimx_1+1:dimx,dimx_1+1:dimx);
      bt=tB1(dimx_1+1:dimx); 
      [t1,t1,Tsi,Tii,noi,conidxi,t1,t1,invTsi]=csd(at,bt,tol*tol);
      beginpi=sum(sum(xindex(1:s-1,:)));
      P1=[P1;TT(beginpi+dims-slambdai1+1:beginpi+dims,:)];
      P2=[P2;TT(beginpi+1:beginpi+dims-slambdai1,:)];
      Tssi=blkdiag(Tssi,Tsi);
      invTssi=blkdiag(invTssi,invTsi);
   end
   Pssi=[P1;P2];
   Tsp=blkdiag(eye(n-ni),Q*Tssi*inv(Pssi));
   Ts=Ts*Tsp;
   invTsp=blkdiag(eye(n-ni),Pssi*invTssi*inv(Q));
   invTs=invTsp*invTs; 
   Tip=blkdiag(eye(m-mi),Toi); 
   Ti=Ti*Tip;
   %at=inv(Ts)*A*Ts;
   %bt=inv(Ts)*B*Ti;
   at=invTs*A*Ts;%%%
   %invTs,B,Ti
   bt=invTs*B*Ti;%%%
   at=at(n-ni+1:n-ni+contr,n-ni+1:n-ni+contr); 
   bt=bt(n-ni+1:n-ni+contr,m-mi+1);
   [t1,t1,tstx,titx,notx,contx]=csd(at,bt,tol*tol);%tol the small the best.
   [t1,t1,tstx,titx,notx,contx,t1,t1,invtstx]=csd(at,bt,tol*tol);%tol the small the best.%%%%%%%%%%%
   ks=[ks contx(1)];
   Tsp=eye(n);
   Tsp(n-ni+1:n-ni+contr,n-ni+1:n-ni+contr)=tstx;
   invTsp=eye(n); invTsp(n-ni+1:n-ni+contr,n-ni+1:n-ni+contr)=invtstx;%%%%
   Ts=Ts*Tsp;
   invTs=invTsp*invTs; 
   Tip=eye(m);
   Tip(m-mi+1,m-mi+1)=titx; 
   Ti=Ti*Tip;
   %at=inv(Ts)*A*Ts;
   %bt=inv(Ts)*B*Ti;
   at=invTs*A*Ts;
   bt=invTs*B*Ti;
   if ni==contr 
      fla=0;
   else 
      a=at(n-ni+contr+1:n,n-ni+contr+1:n);
      b=bt(n-ni+contr+1:n,m-mi+2:m); 
   end
end
%AA=inv(Ts)*A*Ts;
%BB=inv(Ts)*B*Ti; 
AA=invTs*A*Ts;
BB=invTs*B*Ti; 

%Let the elements should be zero be zero.
tt=0*A;
ttt=0*B;
for i=1:length(ks)
   t1=sum(ks(1:i-1));t2=sum(ks(1:i));
   tt(t2,t1+1:t2)=ones(1,t2-t1);
   tt(t1+1:t2-1,t1+2:t2)=eye(ks(i)-1);
   ttt(t2,i)=1;ttt(t1+1:t2,i+1:m)=ones(ks(i),m-i);
end
AA=AA.*tt; 
BB=BB.*ttt; 

%if B is not full column rank.
BB=[BB,zeros(n,mm-m)];
Ti=blkdiag(Ti,eye(mm-m))*U;

err_of_BDCSD_A=norm(AA-invTs*A*Ts);
err_of_BDCSD_B=norm(BB-invTs*BBB*Ti);
