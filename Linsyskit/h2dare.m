function [P,err] = h2dare(A,B,C,D)

%H2DARE  Solution to H2 Discrete-time Algebraic Riccati Equation
%
%        P = h2dare(A,B,C,D)
%
%        returns a positive semi-definite solution, if existent, for
%        the following algebraic Riccati equation for discrete-time
%        H2 optimal control:
%
%            P = A'PA+C'C-(A'PB+C'D)(D'D+B'PB)^{-1}(A'PB+C'D)'
%
%        Note that a positive semi-definite stabilizing solution is
%        existent if and only if the quadruple (A,B,C,D) is left
%        invertible and has no invariant zeros on the unit circle.
%
%        See also DARE, H8DARE and H2CARE.

M=B;
Q=C'*C;
R=D'*D;
N=C'*D;

P=dare(A,M,Q,R,N);
err=inf;

if ~isempty(P)
   err=norm(A'*P*A+C'*C-(A'*P*B+C'*D)*inv(D'*D+B'*P*B)*(A'*P*B+C'*D)'-P);
end

