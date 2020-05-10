function F=dh2state(A,B,C,D,E,tol)

%DH2STATE  Discrete-time H2 Control Using the ATEA Approach
%
%          F = dh2state(A,B,C,D)
%
%          generates a state feedback gain law u(k) = F x(k), which
%          solves H2 optimal control problem for the system:
%
%              x(k+1) = A x(k) + B u(k) + E w(k)
%               h(k)  = C x(k) + D u(k)
%
%          i.e., the H2-norm of the resulting closed-loop transfer
%          matrix from the disturbance, w, to the controlled output,
%          h, is equal to the optimal value, gamma_2^*, which can
%          be pre-calculated using DGM2STAR.
%
%          See also DGM2STAR, H2STATE, DATEA, DH8STATE and DADDPS.

if nargin==4
   tol=1e-8;
end
n=max(size(A)); E=eye(n);
flag=2;
F=datea(A,B,C,D,tol,flag,E);
