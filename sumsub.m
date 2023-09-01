function W = sumsub(X,Y,tol)
% SUMSUB  Orthonormal basis computation for the sum of two subspaces
% W = sumsub(X,Y) returns an orthonormal basis for the sum of two
% subspaces: imX + imY = im[X Y], based on the built-in function 'orth'.
% 
% W = sumsub(X,Y,tol) returns an orthonormal basis for the sum of two
% subspaces: imX + imY = im[X Y], based on the built-in function 'svd',
% where tol is a user-defined threshold to have a certain distance from 
% the machine zero.

% ========================================================================
% Selahattin Burak Sarsilmaz, February 2022 (Matlab R2021a), made use of 
% 1) Function BsumS, Behcet Acikmese,  June 2007.
% Source:
% i) Chapter 3.1.1 of Controlled and Conditioned Invariants in Linear 
% Systems Theory
% ========================================================================

m1 = size(X,1);
m2 = size(Y,1);
if m1 ~= m2
    error('Subspaces do NOT lie in the same vector space.');
end
%% Use default img 
if nargin == 2
    X1 = img(X);
    Y1 = img(Y);
    W = img([X1 Y1]);
end
%% Use svd-based img with tolerance
if nargin == 3
    X1 = img(X,tol);
    Y1 = img(Y,tol);
    W = img([X1 Y1],tol);
end
end
