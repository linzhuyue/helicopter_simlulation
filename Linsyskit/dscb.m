function [AA,BB,CC,DD,Gs,Go,Gi,dims] = dscb(A,B,C,D,tol)

%DSCB  Special Coordinate Basis for Discrete-time Systems
%
%      [As,Bt,Ct,Dt,Gms,Gmo,Gmi,dim] = dscb(A,B,C,D)
%
%      decomposes a discrete-time system characterized by (A,B,C,D)
%      into the standard SCB form with state subspaces x_a being
%      separated into stable, marginally stable and unstable parts
%      (in discrete-time sense), and x_d being decomposed into the
%      form of difference chains.
%
%      Input Parameters:
%
%            x(k+1) = A x(k) + B u(k)
%             y(k)  = C x(k) + D u(k)
%
%      Output Parameters:
%
%            x_t(k+1) = (As+B_0 C_0) x_t(k) + Bt u_t(k)
%             y_t(k)  =       Ct     x_t(k) + Dt u_t(k)
%
%      where x_t = [ x_a^-  x_a^0  x_a^+  x_b  x_c  x_d ]' with
%
%            dim = [ n_a^-, n_a^0, n_a^+, n_b, n_c, n_d ],
%
%      and Gms, Gmo & Gmi = state, output & input transformations.
%
%      See also SCBRAW, SCB and DSSD.

if nargin==4
   tol=1e-8;
end
dc=1;
[AA,BB,CC,DD,Gs,Go,Gi,dims]=scb(A,B,C,D,tol,dc);

%type=2;dc=1;d11_eye=1;
%[AA,BB,CC,DD,Gs,Go,Gi,dims]=zzscbchu(A,B,C,D,tol,type,dc,d11_eye);