function F=dh8state(A,B,C,D,E,gamma,tol)

%DH8STATE  Discrete-time H-infinity Control Using the ATEA Approach
%
%          F = dh8state(A,B,C,D,E,gamma)
%
%          generates a state feedback gain law u(k) = F x(k), which
%          solves H-infinity gamma-suboptimal control problem for
%
%              x(k+1) = A x(k) + B u(k) + E w(k)
%               h(k)  = C x(k) + D u(k)
%
%          i.e., the H-infinity norm of the resulting closed-loop
%          transfer matrix from the disturbance, w, to the controlled
%          output, h, is less than the given 'gamma'. The value of
%          'gamma' has to be chosen larger than the infimum,
%          gamma_\infty^*, which can be determined using DGM8STAR.
%
%          See also DGM8STAR, H8STATE, DATEA, DH2STATE and DADDPS.

if nargin==5
   gamma=1000;tol=1e-8;
elseif nargin==6
   tol=1e-8;
end

flag=8;
F=datea(A,B,C,D,tol,flag,E,gamma);
