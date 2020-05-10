function [Am,Bm,Cm,Dm,Av,Bv,Cv,Dv] = dmpfact(A,B,C,D,tol)

%DMPFACT  Minimum-Phase/All-Pass Factorization of Discrete Systems
%
%         [Am,Bm,Cm,Dm,Av,Bv,Cv,Dv] = dmpfact(A,B,C,D)
%
%         calculates a minimum-phase/all-pass factorization for a
%         detectable system (A,B,C,D) with transfer function matrix
%         G(z), which has no invariant zeros on the unit circle.
%
%         The minimum-phase/all-pass factorization is given as
%
%               G(z) = Gm(z) V(z)
%
%         where
%
%               Gm(z) = Cm (zI - Am)^{-1} Bm + Dm
%
%         is of minimum-phase and left invertible with no infinite
%         zeros, and
%
%               V(z) = Cv (zI - Av)^{-1} Bv + Dv
%
%         is an all-pass factor satisfying V(z) V'(1/z) = I.
%
%         See also DIOFACT and MPFACT.

if nargin==4,
   tol=1e-8;
end;

sys=ss(A,B,C,D,2);
sys=d2c(sys,'tustin');
[A,B,C,D]=ssdata(sys);

[Am,Bm,Cm,Dm,Av,Bv,Cv,Dv]=mpfact(A,B,C,D,tol);

sys=ss(Am,Bm,Cm,Dm);
sys=c2d(sys,2,'tustin');
[Am,Bm,Cm,Dm]=ssdata(sys);
sys=ss(Av,Bv,Cv,Dv);
sys=c2d(sys,2,'tustin');
[Av,Bv,Cv,Dv]=ssdata(sys);