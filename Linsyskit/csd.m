function [AA,BB,Ts,Ti,no,conidx,err_CSD_A,err_CSD_B,invTs]=csd(A,B,tol)

%CSD  Controllability Structural Decomposition for Unsensed Systems
%
%     [At,Bt,Ts,Ti,ucm,Cidx] = csd(A,B)
%
%     returns a controllability structural decomposition for (A,B).
%
%     Input Parameters:
%                 .
%                 x = A x + B u
%
%     Output Parameters:
%
%                           [ Ao  0      0      ...  0     0     ]
%                           | 0   0  I_{k_1-1}  ...  0     0     |
%       At = inv(Ts)*A*Ts = | *   *      *      ...  *     *     |
%                           | :   :      :       .   :     :     |
%                           | 0   0      0      ...  0  I_{k_p-1}|
%                           [ *   *      *      ...  *     *     ]
%
%                           [ 0 ... 0 ]
%                           | 0 ... 0 |
%       Bt = inv(Ts)*B*Ti = | 1 ... 0 |
%                           | :  .  : |
%                           | 0 ... 0 |
%                           [ 0 ... 1 ]
%
%       ucm = uncontrollable modes & Cidx = controllability index
%
%     See also BDCSD, CTRIDX and OSD.

if nargin==2
   tol=1e-8;
end

BBB=B;
mm=size(B,2);
n=size(A,1);

%check if B is full column rank.
[U,t,m]=zzrowup(B',tol);t=t';U=U';
if m==mm
   U=eye(mm);
else
   B=t;
end

if (m==0)|(n==0)
    AA=A;BB=BBB;Ts=eye(n);Ti=eye(mm);no=n;conidx=[];
    err_CSD_A=[];err_CSD_B=[];invTs=eye(n); return;    
end

[AA,BB,Tst,Tit,no,conidx,t1,t1,invTst]=osd(A',B',tol);
tt90=eye(no);
for i=1:length(conidx)
   t1=no+sum(conidx(1:i-1));t2=no+sum(conidx(1:i));
   tt90(t1+1:t2,t1+1:t2)=rot90(eye(conidx(i)));
end
Ti=inv(Tit');
%Ts=inv(Tst')*tt90;
Ts=invTst'*tt90;
invTs=tt90'*Tst';
AA=invTs*A*Ts;
BB=invTs*B*Ti;

%Make the elment should be zero be zero.
tt=0*A; tt(1:no,1:no)=ones(no,no);
ttt=0*B;
if no~=n
   for i=1:length(conidx)
      t1=no+sum(conidx(1:i-1));t2=no+sum(conidx(1:i));
      tt(t2,:)=ones(1,n);
      tt(t1+1:t2-1,t1+2:t2)=eye(conidx(i)-1);
      ttt(t2,i)=1;
   end
end
AA=AA.*tt;
BB=BB.*ttt;

%B is not full column rank
BB=[BB,zeros(n,mm-m)];
Ti=U*blkdiag(Ti,eye(mm-m));%Ti=U*Ti;

err_CSD_A=norm(AA-invTs*A*Ts);
err_CSD_B=norm(BB-invTs*BBB*Ti);
if err_CSD_A>1e-8 | err_CSD_B>1e-8
%   err_CSD_A,err_CSD_B
end
%err_csd=[err_CSD_A,err_CSD_B]
