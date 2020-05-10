function gm2s = gm2s_sc(A,B,E,C1,D1,C2,D2,tol)
%
% gm2s = gm2s_sc(A,B,E,C1,D1,C2,D2,tol)
%
% To calculate the optimal performance value (infimum) for H2 optimization
% for the following continuous-time system
%   .
%   x = A  x + B  u + E  w
%   y = C1 x        + D1 w
%   z = C2 x + D2 u
%
% tol is to determine the result accuracy. Suggest choose tol > 0.001.
%
% This program automatically assume the given system is a singular case!

[n,m] = size(B);
[p,l] = size(D1);

In = eye(n);
Im = eye(m);
Ip = eye(p);

flag = 1;
epsilon = 10;
P = In;
Q = In;
gm2ss = 666;

while flag == 1
    gm2s = gm2ss;
    epsilon = epsilon/10;
    C2t = [C2; epsilon*In; zeros(m,n)];
    D2t = [D2; zeros(n,m); epsilon*Im];
    Et =  [E epsilon*In zeros(n,p)];
    D1t = [D1 zeros(p,n) epsilon*Ip];
    
    P = h2care(A,B,C2t,D2t);
    evP = eig(P);
    if min(evP) < 0
        flag = 0;
    end
    
     Q = h2care(A',C1',Et',D1t');
     evQ = eig(Q);
     if min(evQ) < 0
         flag = 0;
     end

     gm2ss=sqrt(trace(E'*P*E)+trace((A'*P+P*A+C2'*C2)*Q));
  
     if abs(gm2s-gm2ss) < tol
         flag = 0;
     end
end
end

