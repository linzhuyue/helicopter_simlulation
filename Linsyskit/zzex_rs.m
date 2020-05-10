function [A,V,fail]=zzex_rs(A,V,n,l,b1,b2,EPS)
fail=0;
if b1==2,if b2==2,m=l+3;else m=l+2;end;
   x=A(l+1,l+1);y=A(l,l);w=A(l+1,l)*A(l,l+1);p=1.;q=1.;r=1.;
   [A,V]=zzqr_rs(A,V,p,q,r,l,m,n);it=0;it=it+1;
   while it<=30
      z=A(l,l);r=x-z;s=y-z;
      p=(r*s-w)/A(l+1,l)+A(l,l+1);q=A(l+1,l+1)-z-r-s;
      r=A(l+2,l+1);s=abs(p)+abs(q)+abs(r);p=p/s;q=q/s;r=r/s;
      [A,V]=zzqr_rs(A,V,p,q,r,l,m,n);
      if abs(A(m-1,m-2))>EPS*(abs(A(m-1,m-1))+abs(A(m-2,m-2))),it=it+1;
      else A(m-1,m-2)=0.;it=100;
      end;
   end;
   if it<100,fail=1;end;
elseif b2==2,
   x=A(l,l);p=1.;q=1.;r=1.;[A,V]=zzqr_rs(A,V,p,q,r,l,l+2,n);it=0;it=it+1;
   while it<=30
      p=A(l,l)-x;q=A(l+1,l);r=0.;[A,V]=zzqr_rs(A,V,p,q,r,l,l+2,n);
      if abs(A(l+2,l+1))>EPS*(abs(A(l+1,l+1))+abs(A(l+2,l+2))),it=it+1;
      else A(l+2,l+1)=0.;it=100;
      end;
   end;
   if it<100,fail=1;end;
else l1=l+1;q=A(l+1,l+1)-A(l,l);p=A(l,l+1);r=max(abs(p),abs(q));
   if r~=0.
      p=p/r;q=q/r;r=sqrt(p*p+q*q);p=p/r;q=q/r;
      for j=l:n
         s=p*A(l,j)+q*A(l+1,j);A(l+1,j)=p*A(l+1,j)-q*A(l,j);A(l,j)=s;
      end;
      for i=1:l1
         s=p*A(i,l)+q*A(i,l+1);A(i,l+1)=p*A(i,l+1)-q*A(i,l);A(i,l)=s;
      end;
      for i=1:n
         s=p*V(i,l)+q*V(i,l+1);V(i,l+1)=p*V(i,l+1)-q*V(i,l);V(i,l)=s;
      end;
      A(l+1,l)=0.;
   end;
end;

