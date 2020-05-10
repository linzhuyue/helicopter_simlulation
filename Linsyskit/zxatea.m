function Fepsilon = atea(A,B,C,D,option,tol,flag,E,gamma)

%ATEA  Asymptotic Time-scale and Eigenstructure Assignment Design
%
%      F = atea(A,B,C,D[,option])
%
%      produces a state feedback law u = F x using the asymptotic
%      time-scale structure and eigenstructure assignment design
%      method for a continuous-time system characterized by
%              .
%              x = A x + B u,   y = C x + D u
%
%      Users have the 'option' to choose the result either in a
%      numerical or in a symbolic form parameterized by a tuning
%      parameter 'epsilon'. The latter is particularly useful in
%      solving control problems, such as H2 and H-infinity sub-
%      optimal control as well as disturbance decoupling problem.
%      By default or choosing option = 0, the program will ask
%      users to enter a value for 'epsilon' and return a numerical
%      solution. Otherwise, if option = 1, F will be in a symbolic
%      form parameterized by 'epsilon'.
%
%      Note that this function is semi-interactive. Users will be
%      asked to enter desired design eigenstructures during
%      execution. The time-scale is parameterized by the tuning
%      parameter 'epsilon'.
%
%      See also H2STATE, H8STATE, ADDPS and DATEA.

%       flag=0, atea;
%       flag=2, H2;
%       flag=8, H-infinity;
%       flag=100, disturbance decoupling.

clc
disp(' ')
disp('This program will guide your through the step-by-step procedure of the')
disp('Asymptotic Time-scale and Eigenstructure Assignment (ATEA) Design...')

if nargin==4
   option=0;tol=1e-8;flag=0;
elseif nargin==5
   tol=1e-8;flag=0;
elseif nargin==6
   flag=0;
end
Fepsilon=[];

if size(B,2)==0
   disp('...... B is empty......')
   return
end

%STEP ATEA-C.1.
[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=SCB(A,B,C,D,tol);

nan=dims(1);na0=dims(2);nap=dims(3);na=sum(dims(1:3));
nb=dims(4);nc=dims(5);nd=dims(6);
pb=length(lv);mc=length(rv);md=length(qv);
n=size(A,1);
p=size(C,1);
m=size(B,2);
Bd=BB(n-nd+1:n,m0+1:m0+md);

if na0~=0,
   disp('Got invariant zeros on jw axis. I am not programmed to handle such a case...')
   disp(' ')
   return
end
Ass=AA(nan+1:n-nc-nd,nan+1:n-nc-nd);
B0s=BB(nan+1:n-nc-nd,1:m0);
Cd=CC(m0+1:m0+md,n-nd+1:n);
Lsd=AA(nan+1:n-nc-nd,n-nd+1:n)*Cd';
Bs=[B0s Lsd];

Cs=CC(:,nan+1:na+nb);Cs(1:m0,:)=0;
Cs=Go*Cs;
Ds=Go*[eye(m0+md);zeros(p-m0-md,m0+md)];

%STEP ATEA-C.2
[at,bt,t1,T,ks]=ctrbf(Ass,Bs,0*Bs');
ks=sum(ks);
if ks==size(Ass,1)
   at=Ass; bt=Bs; T=eye(ks);
end
P1=eig(at(1:nap+nb-ks,1:nap+nb-ks));
if any(real(P1)>=0)
   disp(' ');
   disp('   System is not stabilizable!');
   return;
end
P2=eig(at(nap+nb-ks+1:nap+nb,nap+nb-ks+1:nap+nb));

Fs=zeros(m0+md,nap+nb);
if length(P2)~=0
   if flag==0
      %clc
      disp(' ');
      disp('Eigenvalue assignment for the slow subsystem: Ass - Bs*Fs ......');
      disp(' ')
      disp(' ')
      disp('The controllable eigenvalues of Ass are given by ...')
      disp(' ')
      disp(P2)
      disp(' ')
      disp('Assign the eigenvalues of the slow subsystem, Ass - Bs*Fs, at ......')
      disp(' ')
      disp('  1). Locations of your choice; or')
      disp('  2). Locations automatically selected by the computer.')
      disp('  ')
      inputno = input('Please select your choice (1 or 2): ');
      if inputno==1
         disp(' ');
         disp(['Please specify ', int2str(length(P2)), ' desired eigenvalues in a row vector.']);
         disp(' ')
         disp(['Note that the desired eigenvalues must have a multiplicity =< ', int2str(rank(bt)), '.']);
         disp(' ')
         disp(' ')
         P3=input('Enter the desired eigenvalues: ');
         Fs=place(at(nap+nb-ks+1:nap+nb,nap+nb-ks+1:nap+nb),bt(nap+nb-ks+1:nap+nb,:),P3);
         Fs=[zeros(m0+md,nap+nb-ks),Fs]*T;
      else
         aass=at(nap+nb-ks+1:nap+nb,nap+nb-ks+1:nap+nb);
         bbss=bt(nap+nb-ks+1:nap+nb,:);
         xss=are(aass,bbss*bbss',0.0001*eye(max(size(aass))));
         Fs=bbss'*xss;  
         Fs=[zeros(m0+md,nap+nb-ks),Fs]*T;
      end
   end
   if flag~=0
      tE=inv(Gs)*E;
      Es=tE(nan+1:nan+nap+nb,:);
   end
   if flag==2
      [Ps,t,Fs,t] = care(Ass,Bs,Cs'*Cs,Ds'*Ds,Cs'*Ds,eye(size(Cs,2)));
      gamma_2_star=sqrt(trace(Es'*Ps*Es))
   end
   if flag==8
      gm8_star=gm8star(Ass,Bs,Cs,Ds,Es,tol)
      disp(' ')
      gamma
      disp(' ')
      while gamma<gm8_star
         gamma = input('Enter the value of gamma, which has to be larger than gamma^*; gamma =  ');
         disp(' ')
      end         
      disp(' ')
      ddt=inv(Ds'*Ds);
      bbt=Es*Es'/gamma/gamma-Bs*ddt*Bs';
      aat=Ass-Bs*ddt*Ds'*Cs;
      cct=Cs'*Cs-Cs'*Ds*ddt*Ds'*Cs;
      Ps=are(aat,-bbt,cct);
      Fs=ddt*(Bs'*Ps+Ds'*Cs);
   end
   if flag==100
      if norm(Es)>tol
         disp(' ')
         disp('The disturbance decoupling problem with state feedback has no solution...')
         disp(' ')
         return;
      else
         disp(' ')
         disp('The disturbance decoupling problem with state feedback has a solution...')
         disp(' ')
         P3=P2-max(real(P2))-1;Fs=place(at(nap+nb-ks+1:nap+nb,nap+nb-ks+1:nap+nb),bt(nap+nb-ks+1:nap+nb,:),P3);
         Fs=[zeros(m0+md,nap+nb-ks),Fs]*T;
      end
   end
end

Fs1=Fs(m0+1:m0+md,:);
Fa1p=Fs(m0+1:m0+md,1:nap);
Fb1 =Fs(m0+1:m0+md,nap+1:nap+nb);
Fa0p=Fs(1:m0,1:nap);
Fb0=Fs(1:m0,nap+1:nap+nb);
Fs0=[Fa0p,Fb0];

%STEP ATEA-C.3
Fc=[];Bc=[];
if nc~=0
   Acc=AA(n-nc-nd+1:n-nd,n-nc-nd+1:n-nd);
   Bc=BB(n-nc-nd+1:n-nd,m-mc+1:m);
   Fc=zeros(mc,nc);
   %clc
   disp(' ');
   disp(' ');
   disp('Eigenvalue Assignment for the slow subsystem: Acc - Bc*Fc: Assign its eigenvalues at');
   disp(' ');
   disp('  1). Locations of your choice; or')
   disp('  2). Locations automatically selected by the computer.')
   disp(' ')
   inputno = input('Enter your option (1 or 2): ');
   disp(' ')
   if inputno==1
      disp(' ');
      disp(['Please specify ', int2str(nc), ' desired eigenvalue(s) in a row vector.']);
      disp(' ');
      disp(['Note that the desired eigenvalues must have a multiplicity =< ', int2str(mc), '.']);
      disp(' ');
      P3=input('Enter the desired eigenvalues: ');
      Fc=place(Acc,Bc,P3);
   else
      Pcc=are(Acc,Bc*Bc',0.0001*eye(max(size(Acc))));
      Fc=Bc'*Pcc;  
   end
end;

%STEP ATEA-C.4
Ft=zeros(md,max(qv));
if md~=0
   %clc
   disp(' ');
   disp('Eigenstructure assignment for fast subsystems, x_{d}, ......');
   disp(' ');
   disp('  1). Specify your own structures; or')
   disp('  2). Let me do it for you.')
   disp(' ')
   inputno = input('Select your option (1 or 2): ');
   if inputno==1
      disp(' ');
      disp('Enter desired eigenvalues for each fast subsystem. The actual closed-loop')
      disp('eigenvalues will be placed at [ the given eigenvalues / epsilon ] ...')
      disp(' ');
      for kk=1:md
         disp(' ');
         disp(['Fast Subsystem No: ', int2str(kk), ', q_',  int2str(kk), ' = ', int2str(qv(kk))])
         disp(' ');
         P3(kk,1:qv(kk))=input(['Enter ', int2str(qv(kk)), ' eigenvalues in row vector: ']);
         tt=poly(P3(kk,1:qv(kk)));
         Ft(kk,1:qv(kk))=tt(2:qv(kk)+1);
      end
   else
      for kk=1:md
         P3=-rand(qv(kk),1);
         P3(1)=-0.5;
         P3=P3-0.5;
         tt=poly(P3);
         Ft(kk,1:qv(kk))=tt(2:qv(kk)+1);
      end
   end
end

if option==1
   syms epsilon
   tFd=sym([]);
   tFa1=sym(zeros(md,nap));
   tFb1=sym(zeros(md,nb));
   tF=sym([]);
   tF0=sym([]);
else
   disp(' ')
   epsilon=input('Enter the value of epsilon: epsilon = ');
   tFd=[];
   tFa1=zeros(md,nap);
   tFb1=zeros(md,nb);
   tF=[];
   tF0=[];
end

for kk=1:md
   for j=1:qv(kk)
      tFi(kk,j)=Ft(kk,qv(kk)-j+1)/epsilon^(qv(kk)-j+1);
   end
end

%STEP ATEA-C.5
for kk=1:md
   if size(Fa1p,2)~=0
      tFa1(kk,:)=Fa1p(kk,:)*tFi(kk,1);
   end
   if size(Fb1,2)~=0
      tFb1(kk,:)= Fb1(kk,:)*tFi(kk,1);
   end
   tFd=blkdiag(tFd,tFi(kk,1:qv(kk)));
end
tFs1=[tFa1 tFb1];

if m0~=0
   tF=[zeros(m0,nan),Fs0,zeros(m0,nc+nd)];
   tF0=CC(1:m0,:);
end
if md~=0
   tF=[tF;zeros(md,nan),tFs1,zeros(md,nc),tFd];
   tF0=[tF0;Bd'*AA(n-nd+1:n,:)];
end
if mc~=0
   tF=[tF;zeros(mc,nan+nap+nb),Fc,zeros(mc,nd)];
   tF0=[tF0;Bc'*AA(n-nc-nd+1:n-nd,1:nan+nap),zeros(mc,nb+nc+nd)];
end

FF=tF+tF0;

Fepsilon=-Gi*(tF+tF0)*inv(Gs);

dig=16;
Fepsilon=vpa(Fepsilon,dig);

%STEP ATEA-C.6
if option==0
   Fepsilon=subs(Fepsilon);
end