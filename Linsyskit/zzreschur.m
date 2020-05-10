function [A,V,Nn,No,Np,ER,EI,TYPE]=zzreschur(A,EPS);
 
%   [A,V,nn,no,np]=zzreschur(A,tol) transforms a square matrix into a real
%	 Schur form, and implements the ordering of the eigenvalues with 
%        the one having the smallest (biggest) real part at the top 
%        (bottom) of the matrix. 
%
%   Input Parameters:
%        A    -  Arbitrary square real matrix to be transformed. 
%        tol  -  Convergence criterio
%
%   Output Parameters:
%	 A    -  Ordered real Schur form matrix
%	 V    -  Orthonormal transformation performed on the input A 
%	 nn   -  Number of eigenvalues with negative real parts
%	 no   -  Number of eigenvalues with real parts equal to zero
%	 np   -  Number of eigenvalues with positive real parts
%

%	 For details on the algorithms, refer to:
%	 G.W. Stewart, "HQR3 and EXCHNG: Fortran Subroutines for Calculationg 
%	 and Ordering the Eigenvalues of a Real Upper Hessenberg Matrix", ACM
% 	 Transactions on Mathematical Software, Vol.2, No.3, 275-280.

%   This program was written by Bo Xiong.

if nargin==1
   EPS=1e-10;
end
[V,A]=schur(A);
n=max(size(A));for i=1:n,TYPE(i)=-1;end;nu=n;
while nu>=1
   l=nu;
   if l>1,
      while abs(A(l,l-1))>EPS*(abs(A(l-1,l-1))+abs(A(l,l))),l=l-1;
         if l==1,break,end;
      end;
   end;
   if l==nu,nl=0;if nu~=1,A(nu,nu-1)=0.;end;TYPE(nu)=0;mu=nu;
      while mu>0,
         if mu==n,mu=nl;nl=0;
	 elseif mu==n-1,
	    if A(mu,mu)<=A(mu+1,mu+1),mu=nl;nl=0;
	    else [A,V,fail]=zzex_rs(A,V,n,mu,1,1,EPS);mu=mu+1;
	    end;
	 elseif abs(A(mu+2,mu+1))<=EPS*(abs(A(mu+1,mu+1))+abs(A(mu+2,mu+2)))
	    A(mu+2,mu+1)=0.;
	    if A(mu,mu)<=A(mu+1,mu+1),mu=nl;nl=0;
	    else [A,V,fail]=zzex_rs(A,V,n,mu,1,1,EPS);mu=mu+1;
	    end;
	 else t1=A(mu,mu);t2=(A(mu+1,mu+1)+A(mu+2,mu+2))/2.;
	    if (t1<=t2)|(abs(t1)+abs(t2)<=EPS),mu=nl;nl=0;
            else [A,V,fail]=zzex_rs(A,V,n,mu,1,2,EPS);
	       if fail==1,TYPE(mu)=-1;TYPE(mu+1)=-1;TYPE(mu+2)=-1;mu=0;l=1;
	       else mu=mu+2;
	       end;
	    end;
         end;
      end;
      nu=l-1;
   else TYPE(nu)=0;TYPE(nu-1)=0;mu=nu;if nu~=2,A(nu-1,nu-2)=0.;end;
      while mu>0,nl=mu-1;
         if mu==n,mu=0;nu=l-1;
         elseif mu==n-1,
            t1=(A(mu-1,mu-1)+A(mu,mu))/2.;t2=A(mu+1,mu+1);
            if (t1<=t2)|(abs(t1)+abs(t2)<=EPS),mu=0;nu=l-1;
            else [A,V,fail]=zzex_rs(A,V,n,nl,2,1,EPS);
               if fail==1,TYPE(nl)=-1;TYPE(nl+1)=-1;TYPE(nl+2)=-1;mu=0;nu=0;
               else mu=mu+1;
               end;
	    end;
	 elseif abs(A(mu+2,mu+1))<=EPS*(abs(A(mu+1,mu+1))+abs(A(mu+2,mu+2)))
	    A(mu+2,mu+1)=0.;
	    t1=(A(mu-1,mu-1)+A(mu,mu))/2.;t2=A(mu+1,mu+1);
	    if (t1<=t2)|(abs(t1)+abs(t2)<=EPS),mu=0;nu=l-1;
	    else [A,V,fail]=zzex_rs(A,V,n,nl,2,1,EPS);
	       if fail==1,TYPE(nl)=-1;TYPE(nl+1)=-1;TYPE(nl+2)=-1;mu=0;nu=0;
	       else mu=mu+1;
	       end;
	    end;
	 else t1=A(mu-1,mu-1)+A(mu,mu);t2=A(mu+1,mu+1)+A(mu+2,mu+2);
	    if (t1<=t2)|(abs(t1)+abs(t2)<=EPS),mu=0;nu=l-1;
            else [A,V,fail]=zzex_rs(A,V,n,nl,2,2,EPS);
	       if fail==1,TYPE(nl)=-1;TYPE(nl+1)=-1;TYPE(nl+2)=-1;TYPE(nl+3)=-1;
	          mu=0;nu=0;
	       else mu=mu+2;
	       end;
	    end;
         end;
      end;
   end;
end;
nu=n;
while nu>=1
   if TYPE(nu)~=-1
      if nu==1,
         ER(nu)=A(nu,nu);EI(nu)=0.;nu=nu-1;
        elseif A(nu,nu-1)==0.,
         ER(nu)=A(nu,nu);EI(nu)=0.;nu=nu-1;
        else 
         [A,V,E1,E2]=zzsp_rs(A,V,n,nu-1);
         if A(nu,nu-1)==0.,ER(nu)=A(nu,nu);EI(nu)=0.;nu=nu-1;
           else 
            ER(nu)=E1;EI(nu-1)=E2;ER(nu-1)=ER(nu);EI(nu)=-EI(nu-1);
            TYPE(nu-1)=1;TYPE(nu)=2;nu=nu-2;
	 end;
      end;
    else nu=nu-1;
   end;
end
for j=1:n-2
    for i=j+2:n,
        A(i,j)=0.;
    end;
end;
Nn=0;No=0;Np=0;
for i=1:n
    if ER(i)<-EPS
       Nn=Nn+1;
       elseif ER(i)>EPS
              Np=Np+1;
       else No=No+1;
    end
end