function [P,flag]=h8dare(A,B,C,D,E,gamma,tol)

%H8DARE  Solution to H-infinity Discrete Algebraic Riccati Equation
%
%        P = h8dare(A,B,C,D,E,gamma)
%
%        returns a positive semi-definite solution, if existent, for
%        the following algebraic Riccati equation for discrete-time
%        H-infinity control:
%
%                           [ B'PA+D'C ]'        [ B'PA+D'C ]
%          P = A'PA + C'C - |          |  G^{-1} |          |
%                           [   E'PA   ]         [   E'PA   ]
%
%        where
%                       [ D'D+B'PB        B'PE     ]
%                   G = |                          |
%                       [   E'PB    E'PE-gamma^2 I ]
%
%        This DARE is related to H-infinity control for the following
%        discrete-time system:
%
%               x(k+1) = A x(k) + B u(k) + E w(k)
%                h(k)  = C x(k) + D u(k)
%
%        Note that a positive semi-definite stabilizing solution is
%        existent if and only if the quadruple (A,B,C,D) is left
%        invertible and has no invariant zeros on the unit circle,
%        as well as gamma > gamma_\infty^*, which can be computed
%        using DGM8STAR.
%
%        See also DARE, H2DARE, H8CARE and DGM8STAR.

if nargin==6
   tol=1e-8;
end

[n,m]=size(B);tm=-inf;a=1;

if isempty(A)
   P=[];
   flag=0;
   return;
end

for k=1:20
   Ft=2*a*rand(m,n)-a;
   t1=min(abs(eig(A+B*Ft)+1));
   if tm<t1
      tm=t1;F=Ft;
   end
end

C2=C;D2=D;
t1=inv(A+B*F+eye(n));
tA=t1*(A+B*F-eye(n));
tB=2*t1*t1*B;
tE=2*t1*t1*E;
tC2=C2+D2*F;
tD2=D2-(C2+D2*F)*t1*B;
tD22=-(C2+D2*F)*t1*E;
tt=[tD2'*tD2,tD2'*tD22;tD22'*tD2,tD22'*tD22-gamma*gamma*eye(size(tD22,2))];
if rank(tt,tol)<size(tt,1)
   P=[];flag=0;
   return;
end
tG1=inv(tt);

BE=[tB'; tE'];
DD=[tD2'; tD22']*tC2;
At=tA-BE'*tG1*DD;
Bt=BE'*tG1*BE;
Ct=tC2'*tC2-DD'*tG1*DD;
Pt=are(At,Bt,Ct);

P=2*t1'*Pt*t1;

flag=(min(real(eig(P)))>-tol);
if flag==0
   P=[];
end