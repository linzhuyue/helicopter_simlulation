function [U,AA,ra]=zzrowup(A,tol,fix)

%rowup [U,AA,ra]=zzrowup(A[,tol])
%     U*A=[ AA ]
%         [  0 ]
%     where AA is full row rank, rank(AA)=ra.
%     and U is unitary matrix.

if nargin==1,
   tol=eps;
end
if nargin<3
   fix=-1;   
end

[U,AA,E]=qr(A); ra=0;
AA=AA*E';
for i=1:min(size(AA)),
   if any(abs(AA(i,:))>tol)
      ra=ra+1;
   end;
end;
if fix~=-1
   ra=fix;
end
AA=AA(1:ra,:);
U=U';