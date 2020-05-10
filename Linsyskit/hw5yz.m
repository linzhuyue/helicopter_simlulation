clc,clear
A=[0,1;0,-1];
B=[0;1];
C1 = [1 0];
D1 = [0 1];
E=[0,0;1,0];
C2=[0,1];
D2=0;
emsilon=0.01;
C2t=[0,1;emsilon,0;0,emsilon;0,0];             
D2t=[0;0;0;emsilon]; 
tol=0.001;
% gm2s = gm2s_sc(A,B,E,C1,D1,C2,D2,tol);
gammaoptimal=0.6185;
% gm8s_sc(A,B,E,C1,D1,C2,D2,gammaoptimal);
P = h2care(A,B,C2t,D2t);                         % right
F = h2state(A,B,C2t,D2t,E);                      % right
F=double(F);
Q = h2care(A',C1',E',D1')';                      % right
K = h2state(A',C1',E',D1',E')';                    % right
K=double(K);
Acmp=A+B*F+K*C1;
Bcmp = -K;
Ccmp = F;
Dcmp=zeros(size(Ccmp,1),size(Bcmp,2));% why??
% closed loop system
Acl = [A+B*Dcmp*C1 B*Ccmp; Bcmp*C1 Acmp];
Bcl = [E+B*Dcmp*D1; Bcmp*D1];
Ccl = [C2+D2*Dcmp*C1 D2*Ccmp];
Dcl = D2*Dcmp*D1;

% Acl = A+B*F;
% Bcl = E;
% Ccl = C2+D2*F;
% Dcl = 0;
w=0.01:0.01:10^4;
Tzw = Tzwo(A,B,E,C1,D1,C2,D2,Acmp,Bcmp,Ccmp,Dcmp,w);
plot(w,Tzw)
xlabel('Frequency(rad/s)');
ylabel('Magnitude)');
set(gca,'XScale','log')
axis([0.1 10^4,-inf,inf])
grid on

[Num,Den] = ss2tf(Acl,Bcl,Ccl,Dcl,size(Ccl,1));
Trans_Func = tf(Num,Den);

figure
bodemag(Trans_Func,{0.01 10000}),grid
alllines = findobj('Type','line');
nlines = length(alllines);
thepoints = cell(nlines,2);
for K = 1:nlines
  thepoints{K,1} = get(alllines(K),'XData');
  thepoints{K,2} = get(alllines(K),'YData');
end
x=thepoints{4,1};y=10.^(thepoints{4,2}./20);
plot(x,y)
grid on
% %-------------------------------------------------------------------------
% clc,clear
% A=[0,1;0,-1];
% B=[0;1];
% C1 = [1 0];
% D1 = [0 1];
% E=[0,0;1,0];
% C2=[0,1];
% D2=0;
% emsilon=0.01;
% C2t=[0,1;emsilon,0;0,emsilon;0,0];             
% D2t=[0;0;0;emsilon]; 
% gamma=2;
% P = h8care(A,B,C2t,D2t,E,gamma)                         % right
% eig(P)
% F = -inv(D2t'*D2t)*(D2t'*C2t+B'*P);
% F=double(F)                                            % double is very important.....
% Q = h8care(A',C1',E',D1',C2',gamma)'                      % right, why the last C2'??
% eig(Q)
% K = -(Q*C1'+E*D1')*inv(D1*D1');                  % right
% K=double(K)                                             % double is very important.....
% Acmp=A+gamma^(-2)*E*E'*P+B*F+inv(eye(2)-gamma^(-2)*Q*P)*K*(C1+gamma^(-2)*D1*E'*P)
% Bcmp=-inv(eye(2)-gamma^(-2)*Q*P)*K
% Ccmp=F
% Dcmp = zeros(size(Ccmp,1),size(Bcmp,2))
% 
% Acl = [A+B*Dcmp*C1 B*Ccmp; Bcmp*C1 Acmp];
% Bcl = [E+B*Dcmp*D1; Bcmp*D1];
% Ccl = [C2+D2*Dcmp*C1 D2*Ccmp];
% Dcl = D2*Dcmp*D1;
% w=0.01:0.01:10^4;
% Tzw = Tzwo(A,B,E,C1,D1,C2,D2,Acmp,Bcmp,Ccmp,Dcmp,w);
% plot(w,Tzw)
% xlabel('Frequency(rad/s)');
% ylabel('Magnitude)');
% set(gca,'XScale','log')
% axis([0.1 10^4,-inf,inf])
% grid on



% [Num,Den] = ss2tf(Acl,Bcl,Ccl,Dcl,size(Ccl,1));
% Trans_Func = tf(Num,Den);
% Trans_Func=Trans_Func/(1-Trans_Func)
% figure
% bodemag(Trans_Func,{0.01 10000}),grid on

% alllines = findobj('Type','line');
% nlines = length(alllines);
% thepoints = cell(nlines,2);
% for K = 1:nlines
%   thepoints{K,1} = get(alllines(K),'XData');
%   thepoints{K,2} = get(alllines(K),'YData');
% end
% x=thepoints{4,1};y=10.^(thepoints{4,2}./20);
% plot(x,y)