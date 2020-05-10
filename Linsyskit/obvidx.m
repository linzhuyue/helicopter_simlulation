function [obidx,no]=obvidx(A,C,tol)

%OBVIDX  Observability Index of Matrix Pair (A,C)
%
%        Oidx = obvidx(A,C)
%
%        returns the observability index for an unforced system.
%
%        Input Parameters:
%            .
%            x = A x,   y = C x
%
%        Output Parameters:
%
%            Oidx = observability index of (A,C)
%
%        See also OSD and CTRIDX.

if nargin==2
   tol=1e-8;
end

[tol,tol,tol,tol,no,obidx]=osd(A,C,tol);
