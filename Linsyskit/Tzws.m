function Tzw = Tzws(A,B,E,C2,D2,F,w)

% Tzw = Tzws(A,B,E,C2,D2,F,w)
%
% Returns the singular values of the closed-loop system comprising of the
% given system
%    .  
%    x = A  x + B  u + E  w
%    y =    x      
%    z = C2 x + D2 u 
% 
% and the state feedback control law
%    
%    u = F x
%
% with frequency range specified in w (rad/s).
%
% See also Tzw_out and ltrloops.

% Programmed by Ben M. Chen - April 17, 2020 at CUHK

Tzw=[]; 
j=sqrt(-1);
D22 = C2*E*0;
for i=1:max(size(w));
    s = j*w(i);
    ts = (C2+D2*F)*inv(s*eye(max(size(A)))-A-B*F)*E+D22;
    Tzw = [Tzw svd(ts)];
end
end
