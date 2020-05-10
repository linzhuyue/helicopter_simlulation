function [P,Q] = gm8sos(A,B,E,C1,D1,C2,D2,gamma,epsilon)
%
% [P,Q] = gm8sos(A,B,E,C1,D1,C2,D2,gamma,epsilon)
%
% To determine whether the given gamma is greater or less than the infimum, 
% i.e., the best achievable performance value for H-infinity optimization  
% for the following continuous-time system:
%   .
%   x = A  x + B  u + E  w
%   y = C1 x        + D1 w
%   z = C2 x + D2 u
%
% If it returns an error message, it likely means the given gamma is less
% than the actual infimum. In such a case, increase the value of gamma &
% try it again.
%
% Perform iteratively to obtain an approximation of the infimum.
% 
% For gamma > the optimal value, this function returns solutions P and Q
% for the corresponding singularly perturbed Riccati equations.
%
% See also gm8star and gm2sos.

% Programmed by Ben M. Chen - April 17, 2020

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
evP = eig(P);

Q = h8care(A',C1',Et',D1t',C2',gamma);
evQ = eig(Q);

if min(evP) < 0 | min(evQ) < 0 | max(eig(P*Q)) > gamma^2
    disp(' ')
    disp('   Given gamma < gm8star  :(  increase gamma and try again...')
    disp(' ')
else
    disp(' ')
    disp('   Given gamma > gm8star  :)  decrease gamma and try again...')
    disp(' ')
end
end

