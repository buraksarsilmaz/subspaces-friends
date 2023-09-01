function W = preimg(A,V,tol)
% PREIMG Orthonormal basis computation for the preimage of a subspace under
% a linear transformation
% W = preimg(A,V) returns an orthonormal basis for the inverse image of V
% under the linear map A via default 'kern'. 
%
% W = preimg(A,V,tol) returns an orthonormal basis for the inverse image of
% V under the linear map A via 'kern', where tol is a user-defined 
% threshold to have a certain distance from the machine zero.

% ========================================================================
% Selahattin Burak Sarsilmaz, February 2022 (Matlab R2021a), made use of 
% 1) Function BVinv, Behcet Acikmese,  June 2007.
% Source:
% i) Chapter 3.1.1 of Controlled and Conditioned Invariants in Linear 
% Systems Theory
% ========================================================================

if isempty(A) || isempty(V) 
    error('The entries of preimg must be nonempty matrices.');
end
m1 = size(A,1);
m2 = size(V,1);
if m1 ~= m2
    error('Dimensions do NOT match.');
end
%% Use default kern 
if nargin == 2
% A^(-1)imV = kerB'A, where B is the basis matrix of kerV'.
B = kern(V');
D = B'*A;
W = kern(D);
end
%% Use svd-based kern with tol
if nargin == 3
B = kern(V',tol);
D = B'*A;
W = kern(D,tol);
end
end
