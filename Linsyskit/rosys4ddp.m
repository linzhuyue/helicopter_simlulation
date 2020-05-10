function [A,B,E,C1,D1,C2,D2,K,flag,N]=rosys4ddp(A,B,E,C1,D1,C2,D2,D22,tol)

%ROSYS4DDP  Irreducible Reduced Order System for Disturbance
%           Decoupling with Static Output Feedback
%
%           [Ar,Br,Er,C1r,D1r,C2r,D2r] = rosys4ddp(A,B,E,C1,D1,C2,D2,D22)
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
%           See also DDPCM.
%

%flag=-1, system is not left invertible, continuous;
%flag=0, DDPCM has not solution, stop;
%flag=1, DDPCM has solution, stop;
%flag=2, system can not be further reduced, stop.

%\Sigma P (A,B,C2,D2);
%\Sigma Q (A,E,C1,D1);

if nargin==8
   tol=1e-10;
end

if ~sum(r_invt(A,B,C2,D2,tol))
   [A,B,E,C1,D1,C2,D2,flag,K,N]=zzrosys4ddp(A,B,E,C1,D1,C2,D2,D22,tol);
   K=K+N;
   return;
end
if ~sum(l_invt(A,E,C1,D1,tol))
   A=A';
   t=B;B=C1';C1=t';
   t=E;E=C2';C2=t';
   t=D1;D1=D2';D2=t';
   D22=D22';
   [A,B,E,C1,D1,C2,D2,flag,K,N]=zzrosys4ddp(A,B,E,C1,D1,C2,D2,D22,tol);
   K=K';N=N';
   K=K+N;
   return;
end

gg=0;flag=-1;N=zeros(size(B,2),size(C1,1));
while flag==-1
   gg=gg+1;
   nt=size(A,1);
   [A,B,E,C1,D1,C2,D2,flag,K1,N1]=zzrosys4ddp(A,B,E,C1,D1,C2,D2,D22,tol);
   if nt==size(A,1)
      flag=2;K=[];
      disp('The system cannot be further reduced, you might need to solve a set of nonlinear')
      disp('polynomial equations to determine if the DDPCM problem is solvable.')
      disp(' ')
      disp('Please refer to the following for detail:')
      disp(' ')
      disp('Linear Systems Theory: A Structural Decomposition Approach')
      disp('Ben M. Chen, Zongli Lin, Yacov Shamash, Birkhauser, 2004')
      disp(' ')
      return;
   end
   if flag==1
      if rem(gg+1,2)
         K1=K1';
      end
      K=N+K1;
      return;
   end
   if flag==0
      K=[];
      return;
   end   
   if rem(gg+1,2)
      N1=N1';
   end
   N=N+N1;
   A=A';
   t=B;B=C1';C1=t';
   t=E;E=C2';C2=t';
   t=D1;D1=D2';D2=t';
   D22=zeros(size(C2,1),size(E,2));
end