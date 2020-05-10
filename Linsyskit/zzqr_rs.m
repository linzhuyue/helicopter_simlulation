function [A,V]=zzqr_rs(A,V,p,q,r,nl,nu,n)
nl2=nl+2;for i=nl2:nu,A(i,i-2)=0.;end;
if nl2~=nu,nl3=nl+3;for i=nl3:nu,A(i,i-3)=0.;end;end;
num1=nu-1;
for k=nl:num1
   last=k==num1;
   x=0.;
   if k~=nl
      p=A(k,k-1);q=A(k+1,k-1);r=0.;if last==0,r=A(k+2,k-1);end;
      x=abs(p)+abs(q)+abs(r);
      if x~=0.
         p=p/x;q=q/x;r=r/x;
      end;
   end;
   if (x~=0.)|(k==nl)
      s=sqrt(p*p+q*q+r*r);if p<0.,s=-s;end;
      if k==nl,if nl~=1,A(k,k-1)=-A(k,k-1);end;
      else A(k,k-1)=-s*x;
      end;
      p=p+s;x=p/s;y=q/s;z=r/s;q=q/p;r=r/p;
      for j=k:n
         p=A(k,j)+q*A(k+1,j);
         if last==0,p=p+r*A(k+2,j); A(k+2,j)=A(k+2,j)-p*z;end;
         A(k+1,j)=A(k+1,j)-p*y;A(k,j)=A(k,j)-p*x;
      end;
      j=min(k+3,nu);
      for i=1:j
         p=x*A(i,k)+y*A(i,k+1);
         if last==0,p=p+z*A(i,k+2);A(i,k+2)=A(i,k+2)-p*r;end;
         A(i,k+1)=A(i,k+1)-p*q;A(i,k)=A(i,k)-p;
      end;
      for i=1:n
         p=x*V(i,k)+y*V(i,k+1);
         if last==0,p=p+z*V(i,k+2);V(i,k+2)=V(i,k+2)-p*r;end;
         V(i,k+1)=V(i,k+1)-p*q;V(i,k)=V(i,k)-p;
      end;
   end;
end;
