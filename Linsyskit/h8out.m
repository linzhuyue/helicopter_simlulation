function [F,K,Acmp,Bcmp,Ccmp,Dcmp,EigCL] = h8out(A,B,E,C1,D1,C2,D2,gamma,epsilon)
%
% [F,K,Acmp,Bcmp,Ccmp,Dcmp,EigCL] = h8out(A,B,E,C1,D1,C2,D2,gamma,epsilon)
%
% Given gamma > the optimal value and a suitable positive scalar epsilon, 
% this function returns an H-infinity gamma suboptimal measurement output 
% feedback control law for the following continuous-time system:
%   .
%   x = A  x + B  u + E  w
%   y = C1 x        + D1 w
%   z = C2 x + D2 u
%
% F is the gain for the corresponding state feedback law u = F x
% K is the observer gain matrix
% 
% The measurement output feedback control law is expressed as
%   .
%   v = Acmp v + Bcmp y
%   u = Ccmp v + Dcmp y
%
% EigCL returns the eigenvalues of the resulting closed-loop system. The
% result is invalid if it is unstable.
%
% See also h8state and h2out.

% Programmed by Ben M. Chen - April 17, 2020 at CUHK

[n,m] = size(B);
[p,l] = size(D1);

In = eye(n);
Im = eye(m);
Ip = eye(p);

C2t = [C2; epsilon*In; zeros(m,n)];
D2t = [D2; zeros(n,m); epsilon*Im];
Et =  [E epsilon*In zeros(n,p)];
D1t = [D1 zeros(p,n) epsilon*Ip];

P = h8care(A,B,C2t,D2t,E,gamma);
Q = h8care(A',C1',Et',D1t',C2',gamma);

F = - inv(D2t'*D2t)*(D2t'*C2t+B'*P);
K = - (Q*C1'+Et*D1t')*inv(D1t*D1t');

Acmp = A+Et*Et'*P/gamma^2+B*F+inv(In-Q*P/gamma^2)*K*(C1+D1t*Et'*P/gamma^2);
Bcmp = -inv(In-Q*P/gamma^2)*K;
Ccmp = F;
Dcmp = Ccmp*Bcmp*0;

% the closed-loop system matrix is for self-testing...
Acl = [A+B*Dcmp*C1 B*Ccmp; Bcmp*C1 Acmp]; 
EigCL = eig(Acl);

end

