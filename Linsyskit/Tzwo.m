function Tzw = Tzwo(A,B,E,C1,D1,C2,D2,Av,Bv,Cv,Dv,w)

% Tzw = Tzwo(A,B,E,C1,D1,C2,D2,Acmp,Bcmp,Ccmp,Dcmp,w)
%
% Returns the singular values of the closed-loop system comprising of the
% given system
%    .  
%    x = A  x + B  u + E  w
%    y = C1 x        + D1 w
%    z = C2 x + D2 u 
% 
% and the measurement output feedback control law
%    .
%    v = Acmp v + Bcmp y
%    u = Ccmp v + Dcmp y
%
% with frequency range specified in w (rad/s).
%
% See also Tzw_state and ltrloops.

% Programmed by Ben M. Chen - April 17, 2020 at CUHK

Tzw=[]; 
j=sqrt(-1);
D22 = C2*E*0;
Cv=-Cv; Dv=-Dv;
for i=1:max(size(w))
    s = j*w(i);
    f1 = (Cv*inv(s*eye(max(size(Av)))-Av)*Bv+Dv)*C1;
    f2 = (Cv*inv(s*eye(max(size(Av)))-Av)*Bv+Dv)*D1;
    tzw = (C2-D2*f1)*inv(s*eye(max(size(A)))-A+B*f1)*(E-B*f2)+(D22-D2*f2);
    Tzw = [Tzw svd(tzw)];
end
end
