function B=sa_act(A,C,tol)

%SA_ACT  Structural Assignment via Actuator Selection
%
%        B = sa_act(A,C)
%
%        For a given unsensed system:
%                  .
%                  x = A x,   y = C x
%
%        the function finds an input matrix B such that the resulting system 
%        characterized by (A,B,C) has the pre-specified desired finite and 
%        infinite zeros.
%
%        Note: Users will be prompted to enter desired structural parameters 
%        after the properties of the given pair (A,C) is evaluated.
%
%        See also SA_SEN.

if nargin==2
   tol=1e-8;
end

B=sa_sen(A',C',tol)';