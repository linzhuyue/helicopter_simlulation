function F=h2state(A,B,C,D,E,option,tol)

%H2STATE  Continuous-time H2 Control Using the ATEA Approach
%
%         F = h2state(A,B,C,D,E[,option])
%
%         generates a state feedback gain law u = F x, which solves
%         H2 suboptimal control problem for the following system:
%              .
%              x = A x + B u + E w
%              h = C x + D u
%
%         i.e., the H2-norm of the resulting closed-loop transfer
%         matrix from the disturbance, w, to the controlled output,
%         h, is minimized. Use GM2STAR to calculate the infimum or
%         the best achievable H2 performance, i.e., gamma_2^*.
%
%         Users have the 'option' to choose the result either in a
%         numerical or in a symbolic form parameterized by a tuning
%         parameter 'epsilon'. By default or choosing option = 0,
%         the program will ask users to enter a value for 'epsilon'
%         and return a numerical solution. Otherwise, if option = 1,
%         F will be in a symbolic form parameterized by 'epsilon'.
%
%         See also GM2STAR, DH2STATE, ATEA, H8STATE and ADDPS.

if nargin==5
   option=0;tol=1e-8;
elseif nargin==6
   tol=1e-8;
end

flag=2;
F=atea(A,B,C,D,option,tol,flag,E);
