function [trfu,trfu1,trfu2]=ss2tf_ds(E,A,B,C,D,tol)

%SS2TF_DS Compute the transfer function of a descriptor system
%
%         [TF,TFs,TFp]=ss2tf_ds(E,A,B,C,D[,tol])
%
%         Input Variables: E, A, B, C, D, and tol (optional)
%
%         Outputs: 
%              TF  = transfer function of descriptor system
%              TFs = strictly proper part of the transfer function
%              TFp = polynomial part of the transfer function


if nargin==5
   tol=1e-10;
end

n=size(A,1);
m=size(B,2);
p=size(C,1);

[PEQ,PAQ,P,Q,n1,n2]=ea_ds(E,A,tol);

N=PEQ(n1+1:n,n1+1:n);A1=PAQ(1:n1,1:n1);
PB=P*B;
CQ=C*Q;
if ~isempty(B)
   B1=PB(1:n1,:);
   B2=PB(n1+1:n,:);
else
   B1=zeros(n,0);B2=zeros(n,0);
end
if ~isempty(C)
   C1=CQ(:,1:n1);
   C2=CQ(:,n1+1:n);
else
   C1=zeros(0,1:n1);C2=zeros(0,n1+1:n);
end

syms s

%trfu1=C1*inv(s*eye(size(A1,1))-A1)*B1;
%trfu2=C2*inv(s*N-eye(size(N,1)))*B2+D

n1=size(A1,1);
sEA1=s*eye(n1)-A1;
dsEA1=det(sEA1);
csEA1=sym(zeros(n1));
if n1>1
   for k=1:n1
      for kk=1:n1
         csEA1(kk,k)=(-1)^(k+kk)*det(sEA1([1:k-1,k+1:n1],[1:kk-1,kk+1:n1]));
      end
   end
else
   csEA1=1;
end

nut=C1*csEA1*B1;
if n1~=0
   trfu1=nut/dsEA1;
else
   trfu1=0;
end

n2=size(N,1);
trfu2=D;
NN=eye(n2);
for kk=1:n2
   trfu2=trfu2-C2*NN*B2*s^(kk-1);
   NN=NN*N;
end

trfu=trfu1+trfu2;
trfu1=vpa(trfu1,16);
trfu2=vpa(trfu2,16);
trfu =vpa(trfu,16);