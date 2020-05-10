function K=ddpcm(A,B,E,C1,D1,C2,D2,D22,tol)

%DDPCM  Disturbance Decoupling with Static Output Feedback
%
%       K = ddpcm(A,B,E,C1,D1,C2,D2,D22)
%
%       computes a solution to the disturbance decoupling problem
%       with a constant (static) measurement output feedback for
%       the following system:
%              .
%              x = A  x + B  u +  E  w
%              y = C1 x        + D1  w
%              h = C2 x + D2 u + D22 w
%
%       if the solution is existent. Otherwise, the program will
%       return an empty matrix for K.
%
%       See also ROSYS4DDP.

if nargin==8
   tol=1e-8;
end

[A,B,E,C1,D1,C2,D2,K,flag,N]=rosys4ddp(A,B,E,C1,D1,C2,D2,D22,tol);
