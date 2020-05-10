function F = addps(A,B,C,D,E,option,tol)

%ADDPS  Solution to Continuous-time Disturbance Decoupling Problem
%
%       F = addps(A,B,C,D,E[,option])
%
%       generates a state feedback gain law u = F x, which solves
%       the almost disturbance decoupling problem for the system:
%              .
%              x = A x + B u + E w
%              h = C x + D u
%
%       i.e., the the resulting closed-loop transfer matrix from the
%       disturbance, w, to the controlled output, h, can be made
%       almost zero. The function will return an empty solution if
%       the problem for the given system is not solvable.
%
%       Users have the 'option' to choose the result either in a
%       numerical or in a symbolic form parameterized by a tuning
%       parameter 'epsilon'. By default or choosing option = 0,
%       the program will ask users to enter a value for 'epsilon'
%       and return a numerical solution. Otherwise, if option = 1,
%       F will be in a symbolic form parameterized by 'epsilon'.
%
%       See also ATEA, H2STATE, H8STATE and DADDPS.

if nargin==5
   option=0;tol=1e-8;
elseif nargin==6
   tol=1e-8;
end

flag=100;
F=atea(A,B,C,D,option,tol,flag,E);