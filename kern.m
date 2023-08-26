function W = kern(X,tol)
% KERN  Orthonormal basis computation for the kernel of a matrix.
% W = kern(X) returns an orthonormal basis for the kernel of X based on
% the built-in function 'null'. 
%
% W = kern(X,tol) returns an orthonormal basis for the kernel of X based on 
% the built-in function 'svd', where tol is a user-defined threshold to
% have a certain distance from the machine zero.

% ========================================================================
% Selahattin Burak Sarsilmaz, February 2022 (Matlab R2021a), made use of 
% 1) Function Bker, Behcet Acikmese,  June 2007.
% ========================================================================

[m,n] = size(X);
%% Use null, which relies on svd
if nargin == 1
     W = null(X);
end
%% Use svd with tolerance
if nargin == 2
[U,S,V] = svd(X);
    if m == 1 || n == 1
       sv = S(1,1);      
    else
       sv = diag(S);   % Singular Values of X
    end
    nz = find(sv>tol);  
    r  = length(nz);   % Rank of X
    if r == n
       W = zeros(n,0);
    else
       W = V(:,r+1:n);
    end
end
end