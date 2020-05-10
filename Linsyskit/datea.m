function F=datea(A,B,C,D,tol,flag,E,gamma)

%DATEA  Eigenstructure Assignment for Discrete-time Systems
%
%       F = datea(A,B,C,D)
%
%       produces a state feedback control law, u(k) = F x(k), using
%       the eigenstructure assignment design method for a discrete-
%       time system characterized by
%
%           x(k+1) = A x(k) + B u(k),   y(k) = C x(k) + D u(k)
%
%       Note that this function is semi-interactive. Users will be
%       asked to enter desired design eigenstructures during
%       execution.
%
%       See also DH2STATE, DH8STATE, DADDPS and ATEA.

%       flag=0, datea;
%       flag=2, H2;
%       flag=8, H-infinity;
%       flag=100, disturbance decoupling.

if nargin==4
   tol=1e-8; flag=0;
elseif nargin==5
   flag=0;
end
F=[];   

%STEP ATEA-D.1.
dc=1;
[Ascb,Bscb,Cscb,Dscb,Gs,Go,Gi,dims,lv,rv,qv,m0]=scb(A,B,C,D,tol,dc);

nan=dims(1);na0=dims(2);nap=dims(3);na=nan+na0+nap;
nb=dims(4);nc=dims(5);nd=dims(6);n=sum(dims);
pb=length(lv);mc=length(rv);md=length(qv);
[p,m]=size(D);

if na0~=0,
   disp(' ')
   disp('(A,B,C,D) has invariant zeros on the unit circle. I am not programmed to handle such a case...')
   disp(' ')
   return
end

Gt=eye(n);
Gt=[Gt(1:nan,:);Gt(na+nb+1:na+nb+nc,:);Gt(nan+1:na+nb,:);Gt(n-nd+1:n,:)];
Gs=Gs*Gt';
Ascb=Gt*Ascb*Gt';
Bscb=Gt*Bscb;
Cscb=Cscb*Gt';

Bd=Bscb(n-nd+1:n,m0+1:m0+md);
Bc=Bscb(nan+1:nan+nc,m0+md+1:m);

Ass=Ascb(nan+nc+1:n,nan+nc+1:n);
Bs=Bscb(nan+nc+1:n,1:m0+md);
Cs=[zeros(m0,nap+nb+nd);Cscb(m0+1:p,nan+nc+1:n)];
Cs=Go*Cs;
Ds=zeros(p,m0+md);Ds(1:m0,1:m0)=eye(m0);
Ds=Go*Ds;

%STEP ATEA-D.2
[at,bt,t1,T,ks]=ctrbf(Ass,Bs,0*Bs');
ks=sum(ks);
if ks==size(Ass,1)
   at=Ass; bt=Bs; T=eye(ks);
end
P1=eig(at(1:nap+nb+nd-ks,1:nap+nb+nd-ks));
if any(abs(P1)>=1),
   disp('');
   disp('   System is not stabilizable!');
   disp('');
   return;
end
P2=eig(at(nap+nb+nd-ks+1:nap+nb+nd,nap+nb+nd-ks+1:nap+nb+nd));

Fs=zeros(m0+md,nap+nb+nd);
if length(P2)~=0,
   if flag==0
      disp(' ');
      disp('Please assign eigenvalues of Ass-Bs*Fs');
      disp(' ')
      disp(' ')
      disp('The controllable eigenvalues of Ass are given by ...')
      disp(' ')
      disp(P2)
      disp(' ')
      disp('Assign the eigenvalues of Ass - Bs * Fs at ......')
      disp(' ')
      disp('  1). Locations of your choice; or')
      disp('  2). Locations automatically selected by the computer.')
      disp('  ')
      inputno = input('Please select your choice (1 or 2): ');
    
      if inputno==1
         disp(' ');
         disp(['Specify  ', int2str(length(P2)), '  desired locations in a row vector.']);
         disp(' ')
         disp('Note that desired locations should be inside the unit circle,');
         disp(['and the desired locations must be of multiplicity =< ', int2str(rank(bt)), '.']);
         disp(' ');
         P3=input('Enter the desired eigenvalues: ');
         Fs=place(at(nap+nb+nd-ks+1:nap+nb+nd,nap+nb+nd-ks+1:nap+nb+nd),bt(nap+nb+nd-ks+1:nap+nb+nd,:),P3);
         Fs=[zeros(m0+md,nap+nb+nd-ks),Fs]*T;
      else
         attt=at(nap+nb+nd-ks+1:nap+nb+nd,nap+nb+nd-ks+1:nap+nb+nd);
         btt=bt(nap+nb+nd-ks+1:nap+nb+nd,:);
         Ptt=dlyap(attt,-btt*btt');
         Ftt=btt'*inv(Ptt);
         Fs=[zeros(m0+md,nap+nb+nd-ks),Ftt]*T;
      end
   end
   if flag~=0
      tE=inv(Gs)*E;
      Es=tE(nan+nc+1:n,:);
   end
   if flag==2
      Ps=h2dare(Ass,Bs,Cs,Ds);
      Fs=inv(Ds'*Ds+Bs'*Ps*Bs)*(Bs'*Ps*Ass+Ds'*Cs);
   end
   if flag==8
      dgm8_star=dgm8star(A,B,C,D,E)
      disp(' ')
      gamma
      disp(' ')
      while gamma<dgm8_star
         gamma = input('Please enter a gamma > gamma^*. gamma = ');
      end
      Ps=h8dare(Ass,Bs,Cs,Ds,Es,gamma);
      t1=inv(Es'*Ps*Es-gamma*gamma*eye(size(Es,2)));
      Fs=inv(Bs'*Ps*Bs+Ds'*Ds-Bs'*Ps*Es*t1*Es'*Ps*Bs)*(Bs'*Ps*Ass+Ds'*Cs-Bs'*Ps*Es*t1*Es'*Ps*Ass);
   end
   if flag==100
      if norm(Es)>tol
         disp(' ')
         disp('The disturbance decoupling problem with state feedback has no solution.')
         disp(' ')
         return;
      else
         disp(' ')
         disp('The disturbance decoupling problem with state feedback has a solution.')
         disp(' ')
         for kk=1:length(P2)
            if abs(P2(kk))>=1
               P2(kk)=0.9/abs(P2(kk))*P2(kk);
            end
         end
         Fs=place(at(nap+nb+nd-ks+1:nap+nb+nd,nap+nb+nd-ks+1:nap+nb+nd),bt(nap+nb+nd-ks+1:nap+nb+nd,:),P2);
         Fs=[zeros(m0+md,nap+nb+nd-ks),Fs]*T;
      end
   end
end

%STEP ATEA-D.3
Fc=[];
if nc~=0
   Acc=Ascb(nan+1:nan+nc,nan+1:nan+nc);
   Bc=Bscb(nan+1:nan+nc,m-mc+1:m);
   Fc=zeros(mc,nc);
   disp(' ');
   disp('* Eigenvalue Assignment for Acc - Bc*Fc, which can be assigned anywhere you like.');
   disp(' ');
   disp(['  Specify ', int2str(nc), ' desired eigenvalue(s) in a self-conjugated row vector.']);
   disp(' ');
   disp('  Note: for stability, the desired eigenvalues should be inside unit circle,');
   disp(['  and each one must be of multiplicity =< ', int2str(mc), '.']);
   disp(' ');
   P3=input('  Enter the desired eigenvalues: ');
   Fc=place(Acc,Bc,P3);
end;

%STEP ATEA-D.4
Eda_m=Bd'*Ascb(n-nd+1:n,1:nan);
Eca_m=Bc'*Ascb(nan+1:nan+nc,1:nan);
Edc=Bd'*Ascb(n-nd+1:n,nan+1:nan+nc);
Eca_p=Bc'*Ascb(nan+1:nan+nc,nan+nc+1:n-nb-nd);
F=zeros(m,n);
F(1:m0,:)=Cscb(1:m0,:);
F(m0+1:m0+md,1:nan+nc)=[Eda_m,Edc];
F(m0+md+1:m,1:nan+nc+nap)=[Eca_m,Fc,Eca_p];
F(1:m0+md,nan+nc+1:n)=F(1:m0+md,nan+nc+1:n)+Fs;
F=-Gi*F*inv(Gs);