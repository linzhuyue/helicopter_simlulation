function [Ai,Bi,Ci,Di,Ao,Bo,Co,Do]=iofact(A,B,C,D,tol)

%IOFACT  Inner-Outer Factorization of Continuous-time Systems
%
%        [Ai,Bi,Ci,Di,Ao,Bo,Co,Do] = iofact(A,B,C,D)
%
%        computes an inner-outer factorization for a stable proper
%        transfer function matrix G(s) with a realization (A,B,C,D),
%        where [B' D'] and [C D] are assumed to be of full rank.
%
%        The inner-outer factorization is given as
%
%                G(s) = Gi(s) Go(s)
%
%        where
%
%                Gi(s) = Ci (sI - Ai)^{-1} Bi + Di
%
%        is an inner, and
%
%                Go(s) = Co (sI - Ao)^{-1} Bo + Do
%
%        is an outer.
%
%        See also MPFACT, GCFACT and DIOFACT.

if nargin==4,
   tol=1e-8;
end

Ai=[];Bi=[];Ci=[];Di=[];Ao=[];Bo=[];Co=[];Do=[];
P=eig(A);
if any(real(P)>=0),
   disp(' ');
   disp('The system is not stable...');
   disp(' ');
   return;
end
[Am,Bm,Cm,Dm,Av,Bv,Cv,Dv]=mpfact(A',C',B',D',tol);
Ai=Av';Bi=Cv';Ci=Bv';Di=Dv';
Ao=Am';Bo=Cm';Co=Bm';Do=Dm';