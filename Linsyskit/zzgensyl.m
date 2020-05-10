function [R,L,err_gensyl]=zzgensyl(A,B,C,D,E,F)

%Gensyl [R,L] = zzgensyl(A,B,C,D,E,F)
% 
%    Solves the Sylvester equations:
%
%               A * R + L * B = C
%               D * R + L * E = F
%

[m,n]=size(C);
Z=[kron(eye(n),A),-kron(B',eye(m));kron(eye(n),D),-kron(E',eye(m))];
colC=[];colF=[];
for i=1:n
   colC=[colC;C(:,i)];
   colF=[colF;F(:,i)];
end;
b=[colC;colF];
x=Z\b;
for i=1:n
   R(:,i)=x((i-1)*m+1:i*m);
   L(:,i)=x(m*n+(i-1)*m+1:m*n+i*m);
end;
L=-L;
err_gensyl=max(norm(A*R+L*B-C),norm(D*R+L*E-F));