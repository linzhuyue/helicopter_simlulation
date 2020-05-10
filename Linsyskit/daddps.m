function F=daddps(A,B,C,D,E,tol)

%DADDPS  Solution to Discrete-time Disturbance Decoupling Problem
%
%        F = daddps(A,B,C,D,E)
%
%        generates a state feedback gain law u(k) = F x(k), which
%        solves the disturbance decoupling problem for the system:
%
%              x(k+1) = A x(k) + B u(k) + E w(k)
%               h(k)  = C x(k) + D u(k)
%
%        i.e., the the resulting closed-loop transfer matrix from the
%        disturbance, w, to the controlled output, h, can be made
%        identically zero. The function will return an empty solution
%        if the problem for the given system is not solvable.
%
%        See also DATEA, DH2STATE, DH8STATE and ADDPS.

if nargin==5
   tol=1e-8;
end
flag=100;
F=datea(A,B,C,D,tol,flag,E);
