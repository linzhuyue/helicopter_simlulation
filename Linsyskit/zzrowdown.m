function [U,AA,ra]=zzrowdown(A,tol,fix)

%rowdown [U,AA,ra]=zzrowdown(A[,tol])
%
%        U*A=[  0 ]
%            [ AA ]
%
%     where AA is full row rank, rank(AA)=ra, and U is unitary matrix.

if nargin==1,
   tol=eps;
end
if nargin<3
   fix=-1;   
end

[U,AA,E]=qr(A); ra=0;
AA=AA*E'; n=size(A,1);
for i=1:min(size(A))
   if any(abs(AA(i,:))>tol)
      ra=ra+1;
   end;
end;
if fix~=-1
   ra=fix;
end
AA=AA(1:ra,:);
tt(1:n-ra,ra+1:n)=eye(n-ra);
tt(n-ra+1:n,1:ra)=eye(ra);
U=tt*U';