function Cout=sa_sen(A,B,tol)

%SA_SEN  Structural Assignment via Sensor Selection
%
%        C = sa_sen(A,B)
%
%        For a given unsensed system:
%                  .
%                  x = A x + B u
%
%        the function finds a measurement output matrix C such that
%        the resulting system characterized by (A,B,C) has the pre-
%        specified desired finite and infinite zeros.
%
%        Note: Users will be prompted to enter desired structural
%        parameters after the properties of the given pair (A,B) is
%        evaluated.
%
%        See also SA_ACT.

if nargin==2
   tol=1e-8;
end

clc

[A1,B1,Ts0,Ti0,no,conidx]=csd(A,B);
[n,m]=size(B);

if no>0
   disp(' ')
   disp(['* The given matrix pair has ' int2str(no) ' uncontrollable (or unobservable) modes,'])
   disp('  which will be included as a part of the invariant zeros of (A,B,C). ')
end
disp(' ')
disp('* The given matrix pair also has an INDEX of')
disp(' ')
disp(conidx)

disp('* Assignment of desired infinite zero structure, i.e., Morse List I4 ...')

ff=0;
while ff==0
   disp(' ')
   disp('  I4 List can be assigned subject to the following constrains:')
   disp(' ')
   disp(['   1) the number of items in I4 should be = ', int2str(m), ', and'])
   disp('   2) each entry of I4 must be > 0 but =< the entry contained in the INDEX.')
   disp(' ')   
   I4=input('  Enter a desired I4 = ');
   I4=sort(round(I4));
   ff=1;
   
   if min(size(I4))~=1 | max(size(I4))~=m | min(conidx-I4)<0 | min(I4)<=0
      disp(' ')
      disp('  > > > The I4 List entered is invalid. Please re-enter a valid one. < < <')
      ff=0;
   end
end

tcd=conidx;
tt=eye(n);
P2(1:no,1:n)=tt(1:no,:);
nd=sum(I4);
na=n-nd;
C4=zeros(m,n);
ttd=zeros(0,n);
td=0;
ttd=zeros(0,n);
for k=1:m
   nk_1=no+sum(conidx(1:k-1));
   nk=no+sum(conidx(1:k));
   if tcd(k)~=inf
      td=td+1;
      P2=[P2;tt(nk_1+1:nk_1+conidx(k)-I4(td),:)];
      ttd=[ttd;tt(nk_1+conidx(k)-I4(td)+1:nk,:)];
      C4(td,nk_1+conidx(k)-I4(td)+1)=1;
   end
end
P2=[P2;ttd];
P2=P2';
A4=P2'*A1*P2; B4=P2'*B1; C4=C4*P2;

if na~=0
   Aa=A4(1:na,1:na);
   Lad=A4(1:na,n-nd+1:n)*C4(:,n-nd+1:n)';
   Aa_ast=A4(no+1:na,no+1:na);
   Lad_ast=Lad(no+1:na,:);
   Ao=Aa(1:no,1:no);
   
   disp(' ')
   disp('* Assignment of desired invariant zeros ... ')
   disp(' ')
   disp(['  You are free to assign ', int2str(na-no),' invariant zero(s).']);
   
   if na-no>0
      disp(' ')
      disp(['  You need to enter ', int2str(na-no),' desired invariant zero(s) in a self-conjugated set']);
      disp(['  with each zero having a multiplicity less than or equal to ', int2str(min([m,size(Aa_ast,1),rank(Lad_ast)])), '.'])
      disp(' ')
      Pc2=input('  Enter the desired invariant zero(s) = ');
      Kc=[rand(size(C4,1),no),place(Aa_ast,Lad_ast,Pc2)];
   else
      Kc=rand(m,na);
   end
   C4(1:m,1:na)=Kc;
end

Cout=C4*inv(Ts0*P2);