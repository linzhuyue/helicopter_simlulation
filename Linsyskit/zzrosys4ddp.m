function [A_R,B_R,E_R,C1_R,D1_R,C2_R,D2_R,flag,K,N]=zzrosys4ddp(A,B,E,C1,D1,C2,D2,D22,tol)

%ZZROSYS4DDP  Reduced Order System for Disturbance Decoupling with Static Output Feedback
%
%           [Ar,Br,Er,C1r,D1r,C2r,D2r,flag,K,N] = zzrosys4ddp(A,B,E,C1,D1,C2,D2,D22)
%
%           generates an irreducible reduced order system from the
%           original system:
%                 .
%                 x = A  x + B  u +  E  w
%                 y = C1 x        + D1  w
%                 h = C2 x + D2 u + D22 w
%
%           The reduced order system is characterized by
%                 .
%                 x_r = Ar  x_r + Br  u_r +  Er  w_r
%                 y_r = C1r x_r           + D1r  w_r
%                 h_r = C2r x_r + D2r u_r
%
%           which can be used to solve the static output disturbance
%           decoupling problem for the original system through some
%           numerical computation package such as QEPCAD.
%
%           The program will return an empty result if the problem
%           for the original system is not solvable.
%
%           flag=0, DDPCM has not solution
%           flag=1, DDPCM has solution K
%           flag=-1, the system (A,B,C2,D2) is not left invertible.
%
%           See also DDPCM, ROSYS4DDP

%System P (A,B,C2,D2);
%System Q (A,E,C1,D1);

if nargin==8
   tol=1e-8;
end

[ell,n]=size(C1);
m=size(B,2);
[p,q]=size(D22);

V_star_p=v_star(A,B,C2,D2,tol);
S_star_q=s_star(A,E,C1,D1,tol);

%STEP DDPCM-R.O.S.1
X=zznulltol(V_star_p',tol)';
Y=S_star_q;
t1=[X*B;D2];
t2=[C1*Y,D1];
t3=[X*A*Y,X*E;C2*Y,D22];
t=-pinv(t1)*t3*pinv(t2);
flag=(norm(t1*t*t2+t3)<tol);%Necessary condition in Theorem 11.3.2
%zt=rand(size(t));N=t+zt+pinv(t1)*t1*zt*t2*pinv(t2)
N=t;
%N=-pinv(t1'*t1)*t1'*t3*t2'*pinv(t2*t2');


K=[];
if flag==0
   disp(' ')
   disp('DDPCM has no solution...')
   disp(' ')
   A_R=[];B_R=[];E_R=[];C1_R=[];D1_R=[];C2_R=[];D2_R=[];N=[];
   return;
end

%STEP DDPCM-R.O.S.2
[G_m,D10,nr10]=zzrowup(D1,tol);
G_m=G_m';
C1m=G_m'*C1;
D1m=G_m'*D1;
A_N=A+B*N*C1;
E_N=E+B*N*D1;
C_N=C2+D2*N*C1;

%STEP DDPCM-R.O.S.3
dc=0;
[As_N,Bs_N,Cs_N,Ds_N,Gs,Go,Gi,dims_N,lv_N,rv_N,qv_N,m0_N]=scb(A_N,B,C_N,D2,tol);
%type=2;dc=0;d11_eye=1;
%[As_N,Bs_N,Cs_N,Ds_N,Gs,Go,Gi,dims_N,lv_N,rv_N,qv_N,m0_N]=zzscbchu(A_N,B,C_N,D2,tol,type,dc,d11_eye);

n_N=size(A_N,1); na_N=sum(dims_N(1:3)); nb_N=dims_N(4);
nc_N=dims_N(5);  nd_N=dims_N(6); md_N=length(qv_N);
Gt=eye(n_N); Gt=[Gt(:,na_N+nb_N+1:na_N+nb_N+nc_N),Gt(:,1:na_N+nb_N),Gt(:,n_N-nd_N+1:n_N)];
As_N=Gt'*As_N*Gt;
Bs_N=Gt'*Bs_N;
Cs_N=Cs_N*Gt;
Gs=Gs*Gt;
Bd_N=Bs_N(n_N-nd_N+1:n_N,m0_N+1:m0_N+md_N);

C1=inv(G_m)*C1*Gs;
D1=inv(G_m)*D1;
E_abcd=inv(Gs)*E_N;

%STEP DDPCM-R.O.S.4
Aaa=As_N(nc_N+1:nc_N+na_N,nc_N+1:nc_N+na_N);
Ea=E_abcd(nc_N+1:nc_N+na_N,:);

[A_4,E_4,t,G_a,kt]=ctrbf(Aaa,Ea,0*Ea');
n_ac=sum(kt);
nt=size(Aaa,1);t1=eye(nt);t1=[t1(nt-n_ac+1:nt,:);t1(1:nt-n_ac,:)];
Gt=blkdiag(eye(nc_N),G_a'*t1');
%Gt,tt=inv(G_a'*t1')*Aaa*(t1*G_a)-A_4

%STEP DDPCM-R.O.S.5
nt=nc_N+na_N;
A_R=inv(Gt)*As_N(1:nt,1:nt)*Gt;
B_R=inv(Gt)*Bs_N(1:nt,:);
E_R=inv(Gt)*E_abcd(1:nt,:);
C1_R=C1(:,1:nt)*Gt;
D1_R=D1;
C2_R=Cs_N(:,1:nt)*Gt;
BdEd_R=As_N(n_N-nd_N+1:n_N,1:nt)*Gt;

nt=nc_N+n_ac;
A_R=A_R(1:nt,1:nt);
B_R=B_R(1:nt,:);
E_R=E_R(1:nt,:);
C1_R=C1_R(:,1:nt);
C2_R=C2_R(:,1:nt);
BdEd_R=BdEd_R(:,1:nt);
Cc11a=C1_R(nr10+1:size(C1_R,1),nc_N+1:nc_N+n_ac);

A_R=A_R+B_R(:,1:m0_N)*C2_R(1:m0_N,:);
B_R=B_R*inv(Gi);
C1_R=G_m*C1_R;
D1_R=G_m*D1_R;
C2_R=[C2_R(1:m0_N,:);Bd_N'*BdEd_R];
D2_R=zeros(size(C2_R,1),size(B_R,2)); D2_R(1:m0_N+md_N,1:m0_N+md_N)=eye(m0_N+md_N);
D2_R=D2_R*inv(Gi);
C0aEda=C2_R(:,nc_N+1:nc_N+n_ac);

p_R=size(C0aEda,1);
ell_h=size(Cc11a,1);

rv=r_invt(A,B,C2,D2,tol);
if sum(rv)==0   
   t1=zznulltol(Cc11a',tol);
   t2=zznulltol(C0aEda',tol);
   flag=ssorder(t1,t2);
   flag=(flag<1);
   
   %flag=(norm(C0aEda+K0d*Cc11a)<tol);
   if flag==0
      disp(' ')
      disp(':-(  ~ DDPCM has no solution...')
      disp(' ')
      N=[];
   elseif flag==1
      disp(' ')
      disp(':-)  ~ DDPCM has a solution...')
      disp(' ')
      if isempty(Cc11a)
         K0d=rand(size(Cc11a,1));%K0d=eye(size(Cc11a,1));
      else
         K0d=-C0aEda*pinv(Cc11a);
      end
      K=Gi*[zeros(p_R,ell-ell_h),K0d;rand(m-p_R,ell)]*inv(G_m);
   end
else
   disp('(A,B,C2,D2) is not left invertible...')
   disp(' ')
   flag=-1;
end