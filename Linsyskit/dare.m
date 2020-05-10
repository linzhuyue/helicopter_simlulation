function [P,err] = dare(A,M,Q,R,N)

%DARE  Solution to Discrete-time Algebraic Riccati Equation
%
%      [P,err] = dare(A,M,Q,R,N)
%
%      returns a positive semi-definite solution, if existent, for
%      the following general discrete algebraic Riccati equation:
%
%        P = A'PA - (A'PM + N)(R + M'PM)^{-1}(M'PA + N') + Q
%
%      using non-iterative method reported in Chen's work, Robust
%      and H-infinity Control, Springer, London, 2000.
%
%      err is the solution error defined as:
%
%        err = | A'PA-(A'PM+N)(R+M'PM)^{-1}(M'PA+N')+Q-P |
%
%      which can be used to verify the accuracy of the solution, P,
%      to the DARE.
%
%      See also H2DARE and H8DARE.

P=[];err=[];
if isempty(A)|isempty(M)|isempty(Q)|isempty(R)|isempty(N)
   warning('A,M,Q,R,N cannot be empty matrix')
   return;
end

n=size(A,1);
A_1=inv(A+eye(n));
A_2=inv(A'+eye(n));

F=A_1*(A-eye(n));
G=2*A_1*A_1*M;
W=R+M'*A_2*Q*A_1*M-N'*A_1*M-M'*A_2*N;
H=-Q*A_1*M+N;

if rank(W)<size(W,1)
   P=[];err=inf;
   return;
end   
   
invW=inv(W);
At=F-G*invW*H';
Bt=G*invW*G';
Ct=Q-H*invW*H';

tP=are(At,Bt,Ct);
P=2*A_2*tP*A_1;

err=norm(A'*P*A-(A'*P*M+N)*inv(R+M'*P*M)*(M'*P*A+N')+Q-P);