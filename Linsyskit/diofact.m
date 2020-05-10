function [Ai,Bi,Ci,Di,Ao,Bo,Co,Do] = diofact(A,B,C,D,tol)

%DIOFACT  Inner-Outer Factorization of Discrete-time Systems
%
%         [Ai,Bi,Ci,Di,Ao,Bo,Co,Do] = diofact(A,B,C,D)
%
%         computes an inner-outer factorization for a stable proper
%         transfer function G(z) with a realization (A,B,C,D).
%
%         The inner-outer factorization is given as
%
%                G(z) = Gi(z) Go(z)
%
%         where
%
%                Gi(z) = Ci (zI - Ai)^{-1} Bi + Di
%
%         is an inner, and
%
%                Go(z) = Co (zI - Ao)^{-1} Bo + Do
%
%        is an outer.
%
%        See also DMPFACT and DIOFACT.

if nargin==4
   tol=1e-8;
end

sys=ss(A,B,C,D,2);
sys=d2c(sys,'tustin');
[A,B,C,D]=ssdata(sys);

[Ai,Bi,Ci,Di,Ao,Bo,Co,Do]=iofact(A,B,C,D,tol);

sys=ss(Ai,Bi,Ci,Di);
sys=c2d(sys,2,'tustin');
[Ai,Bi,Ci,Di]=ssdata(sys);
sys=ss(Ao,Bo,Co,Do);
sys=c2d(sys,2,'tustin');
[Ao,Bo,Co,Do]=ssdata(sys);