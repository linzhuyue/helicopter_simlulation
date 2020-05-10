function [AA,BB,CC,DD,Gs,Go,Gi,dim,m0]=scbraw(A,B,C,D,tol)

%SCBRAW  Special Coordinate Basis in a Raw Form
%
%        [As,Bt,Ct,Dt,Gms,Gmo,Gmi,dim] = scbraw(A,B,C,D)
%
%        decomposes a proper system characterized by (A,B,C,D) into
%        a raw SCB form without separating state subspace x_d into
%        integration chains.
%
%        Input Parameters:
%              .
%              x = A x + B u,   y = C x + D u
%
%        Output Parameters:
%              .
%              x_t = At x_t + Bt u_t,   y_t = Ct x_t + Dt u_t
%
%        where x_t = [ x_a  x_b  x_c  x_d ]' with dimensions of
%
%              dim = [ n_a, n_b, n_c, n_d ], respectively,
%
%        and Gms, Gmo & Gmi = state, output & input transformations.
%
%        See also SCB and DSCB.

if nargin==4
   tol=1e-8;
end
type=0;dc=0;d11_eye=1;
[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);
dim=sum(dims(1:3));
dim(2:4)=dims(4:6);