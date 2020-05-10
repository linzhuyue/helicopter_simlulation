function F=h8state(A,B,C,D,E,gamma,option,tol)

%H8STATE  Continuous-time H-infinity Control Using the ATEA Approach
%
%         F = h8state(A,B,C,D,E,gamma[,option])
%
%         generates a state feedback gain law u = F x, which solves
%         H-infinity gamma-suboptimal control problem for the system:
%              .
%              x = A x + B u + E w
%              h = C x + D u
%
%         i.e., the H-infinity norm of the resulting closed-loop
%         transfer matrix from the disturbance, w, to the controlled
%         output, h, is less than the given 'gamma'. The value of
%         'gamma' has to be chosen larger than the infimum,
%         gamma_\infty^*, which can be pre-determined using GM8STAR.
%
%         Users have the 'option' to choose the result either in a
%         numerical or in a symbolic form parameterized by a tuning
%         parameter 'epsilon'. By default or choosing option = 0,
%         the program will ask users to enter a value for 'epsilon'
%         and return a numerical solution. Otherwise, if option = 1,
%         F will be in a symbolic form parameterized by 'epsilon'.
%
%         See also GM8STAR, DH8STATE, ATEA, H2STATE and ADDPS.

if nargin==5
   gamma=100;option=0;tol=1e-8;
elseif nargin==6
   option=0;tol=1e-8;
elseif nargin==7
   tol=1e-8;
end

flag=8;
F=atea(A,B,C,D,option,tol,flag,E,gamma);
