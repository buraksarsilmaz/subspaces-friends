function X = stabsub(A,tol)
% STABSUB Orthonormal basis computation for the stable subspace of a square
% real matrix
% W = stabsub(A) returns an orthonormal basis for the stable subspace 
% of A.
%
% W = stabsub(A,tol) returns an orthonormal basis for the stable subspace 
% of A, where tol is a user-defined threshold to have a certain distance
% from the machine zero. 

% ========================================================================
% Selahattin Burak Sarsilmaz, August 2023 (Matlab R2023a), made use of 
% 1) Function Bstabsub, Behcet Acikmese,  June 2007. 
% Source: 
% i) Theorem 2. 1 and 2.14 of Control theory for linear systems, 
% Trentelman, Stoorvogel, and Hautus.
% ========================================================================

n1 = size(A,1);
n2 = size(A,2);
if n1 ~= n2 || isempty(A)
    error('The first input must be a nonempty square matrix.');
end

v = eig(A);  
i = real(v) < 0; % indices of stable eigenvalues
vs = v(i);       % stable eigenvalues

if length(vs) == n1
    p = zeros(n1,n1);
else 

I = eye(n1);
p = I;
for k = 1 : length(vs)
    p = p*(A-vs(k)*I);
end

    p = real(p);        %since A is real 
end

if nargin < 2
    X = kern(p);
else
    X = kern(p,tol);
end
   
end

