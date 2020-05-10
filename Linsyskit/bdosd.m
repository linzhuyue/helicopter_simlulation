function [AA,CC,Ts,To,ks,err_of_BDOSD_A,err_of_BDOSD_C,invTs,invTo] = bdosd(A,C,tol)
 
%BDOSD  Block Diagonal Observable Structural Decomposition
%
%       [At,Ct,Ts,To,ks] = bdosd(A,C)
%
%       transforms an observable pair (A,C) into the block diagonal
%       observable structural decomposition form.
%
%       Input Parameters:
%              .
%              x = A x,   y = C x
%
%       Output Parameters:
%
%                               [ A_1   0    0  ...   0  ]
%                               |  0   A_2   0  ...   0  |
%           At = inv(Ts)*A*Ts = |  0    0   A_3 ...   0  |
%                               |  :    :    :   .    :  |
%                               [  0    0    0  ...  A_k ]
%
%                               [ C_1   0    0  ...   0  ]
%                               |  *   C_2   0  ...   0  |
%           Ct = inv(To)*C*Ts = |  *    *   C_3 ...   0  |
%                               |  :    :    :   .    :  |
%                               |  *    *    *  ...  C_k |
%                               [  *    *    *  ...   *  ]
%       where
%                               [ *  1  0 ... 0 ]
%                               | *  0  1 ... 0 |
%                         A_i = | :  :  :  .  : |
%                               | *  0  0 ... 1 |
%                               [ *  0  0 ... 0 ]
%
%                         C_i = [ 1  0  0 ... 0 ]
%
%           ks: contains the sizes of blocks A_1, ..., A_k.
%
%       See also OSD and BDCSD.
 
if nargin==2
   tol=1e-8;
end

p=size(C,1);
n=size(A,1);
if (p==0)|(n==0)
   AA=A;CC=C;Ts=eye(n);To=eye(p);ks=0;
   err_of_BDOSD_A=0;err_of_BDOSD_C=0;invTs=eye(n);invTo=eye(p);
   return;
end

[AA,CC,ts,to,ks,t1,t1,flag]=bdcsd(A',C',tol);
if flag==0
   AA=A;CC=C;
   Ts=eye(size(A,1));To=eye(size(C,2));ks=[];
   err_of_BDOSD_A=0;err_of_BDOSD_B=0;invTs=Ts;invTo=To;
   return;
end

Ts=inv(ts'); invTs=ts';
To=inv(to'); invTo=to'; 
 
tt=[]; 
for i=1:length(ks) 
   t1=sum(ks(1:i-1));t2=sum(ks(1:i)); 
   tt(t1+1:t2,t1+1:t2)=rot90(eye(ks(i))); 
end 
Ts=Ts*tt;
invTs=tt'*invTs; 
AA=invTs*A*Ts; 
CC=invTo*C*Ts; 

% Make the elment should be zero be zero. 
tt=zeros(size(A)); 
ttt=zeros(size(C)); 
for i=1:length(ks) 
   t1=sum(ks(1:i-1));t2=sum(ks(1:i)); 
   tt(t1+1:t2,t1+1)=ones(t2-t1,1); 
   tt(t1+1:t2-1,t1+2:t2)=eye(ks(i)-1); 
   ttt(i,t1+1)=1; ttt(i+1:p,t1+1:t2)=ones(p-i,ks(i)); 
end 
AA=AA.*tt; 
CC=CC.*ttt; 
 
err_of_BDOSD_A=norm(AA-invTs*A*Ts);
err_of_BDOSD_C=norm(CC-invTo*C*Ts);
%if err_of_BDOSD_A>1e-8
%   err_of_BDOSD_A, err_of_BDOSD_C
%end
