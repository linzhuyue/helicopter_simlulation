function [Lt,La,E]=ltrloops(A,B,C,D,F,Ac,Bc,Cc,Dc,w)
%
% [Lt,La,E]=ltrloops(A,B,C,D,F,Acmp,Bcmp,Ccmp,Dcmp,w) 
%
% generates the following frequency responses for Loop Transfer Recovery:
%
%  Lt - Singular value of targetr loop
%  La - Singular value of achieved loop with output feedback controller
%  E  - Singular value of recovery error
%
% for a given system
%  .
%  x = A x + B u
%  y = C x + D u
%
% where F is the target state feedback gain (i.e., u = - F x) & the output 
% the output feedback controller is characterized by
%  .
%  v =   Acmp v + Bcmp y
%  u = - Ccmp v - Dcmp y
%
%  w contains the frequencies in rad/s.
%
% See also Tzw_state and Tzw_out.
 
% Programmed by Ben M. Chen 
% School of Electrical Engineering and Computer Science
% Washington State University 
% Pullman, WA 99164-2752.
%   
% July 9, 1991
%
% Modified by Ben M. Chen on April 17, 2020 at CUHK
 
[n,m]=size(A);
In=eye(n);
[nc,x]=size(Ac);
Inc=eye(nc);
Lt=[];La=[];E=[];
[wm,wn]=size(w);
for i=1:max(wm,wn)   
    s=sqrt(-1)*w(i);
    Lti=F*inv(s*In-A)*B;
    Lai=(Cc*inv(s*Inc-Ac)*Bc+Dc)*(C*inv(s*In-A)*B+D);
    Ei =Lti-Lai; 
    Lt =[Lt svd(Lti)];
    La =[La svd(Lai)];
    E  =[E  svd(Ei)];
end
end