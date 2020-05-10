function [AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0,err_scbbest,for_gm8_Gmor,tol]=SCB(A,B,C,D,tol,dc)

%SCB  Special Coordinate Basis for Continuous-time Systems
%
%     [As,Bt,Ct,Dt,Gms,Gmo,Gmi,dim] = scb(A,B,C,D)
%
%     decomposes a continuous-time system characterized by (A,B,C,D)
%     into the standard SCB form with state subspaces x_a being
%     separated into stable, marginally stable and unstable parts
%     (in continuous-time sense), and x_d being decomposed into the
%     form of integration chains.
%
%     Input Parameters:
%           .
%           x = A x + B u,   y = C x + D u
%
%     Output Parameters:
%           .
%           x_t = (As + B_0 C_0) x_t + Bt u_t,   y_t = Ct x_t + Dt u_t
%
%     where x_t = [ x_a^-  x_a^0  x_a^+  x_b  x_c  x_d ]' with
%
%           dim = [ n_a^-, n_a^0, n_a^+, n_b, n_c, n_d ],
%
%     and Gms, Gmo & Gmi = state, output & input transformations.
%
%     See also SCBRAW, DSCB and SSD.

%     dc : dc=0, for continuous-time system; (default)
%          dc=1, for discrete-time system.

%     Note: make a choice between chen's and chu's algorithm. The smaller error one is selected.


if nargin==4,
   tol=1e-8; dc=0;   
elseif nargin==5
   dc=0;   
end

type=2;d11_eye=1;
[AA1,BB1,CC1,DD1,Gs1,Go1,Gi1,dims1,lv1,rv1,qv1,m01,err_scb,for_gm8_Gmor1]=zzscb(A,B,C,D,tol,dc);
[AA,BB,CC,DD,Gs,Go,Gi,dims,lv,rv,qv,m0,err_scbbest,for_gm8_Gmor]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);
if max(err_scb)<max(err_scbbest)
   AA=AA1;
   BB=BB1;
   CC=CC1;
   DD=DD1;
   Gs=Gs1;
   Go=Go1;
   Gi=Gi1;
   lv=lv1;
   rv=rv1;
   qv=qv1;
   m0=m01;
   dims=dims1;
   err_scbbest=err_scb;
   for_gm8_Gmor=for_gm8_Gmor1;
end