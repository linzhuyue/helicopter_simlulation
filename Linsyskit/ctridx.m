function [Cidx,no]=ctridx(A,B,tol)

%CTRIDX  Controllability Index of Matrix Pair (A,B)
%
%        Cidx = ctridx(A,B)
%
%        returns the controllability index for an unsensed system.
%
%        Input Parameters:
%            .
%            x = A x + B u
%
%        Output Parameters:
%
%            Cidx = controllability index of (A,B)
%
%        See also CSD and OBVIDX.

if nargin==2
   tol=1e-8;
end

[tol,tol,tol,tol,no,Cidx]=csd(A,B,tol);
