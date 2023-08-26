function W = isalci(K,A,B,tol)
% ISALCI Orthonormal basis computation for the largest (A,B)-controlled
% invariant subspace in the image of K. 
% W = isalci(K,A,B) returns an orthonormal basis for the the largest 
% (A,B)-controlled invariant subspace in the image of K through the
% invariant subspace algorithm. 
% 
% W = isalci(K,A,B,tol) returns an orthonormal basis for the the largest 
% (A,B)-controlled invariant subspace in the image of K through the 
% invariant subspace algorithm, where tol is a user-defined threshold to
% have a certain distance from the machine zero. 

% ========================================================================
% Selahattin Burak Sarsilmaz, February 2022 (Matlab R2021a), made use of 
% 1) Function BciL, Behcet Acikmese,  June 2007. 
% Sources: 
% i) Theorem 4.10 of Control theory for linear systems, 
% Trentelman, Stoorvogel, and Hautus.
% ii) Algorithm 4.1.2 of Controlled and conditioned invariants in linear
% systems theory, Basile and Marro.  
% ========================================================================

if isempty(A) || isempty(B) || isempty(K)
    error('The entries of isalci must be nonempty matrices.');
end
[nx1,nx2] = size(A);
nx3 = size(B,1);
nv1 = size(K,1);
if nx1 ~= nx2 || nx2 ~= nx3 || nv1 ~= nx2
    error('Dimensions do NOT match.')
end

if nargin < 4
    V = img(K);
else
    V = img(K,tol);
end
n = size(V,2);
flag = 0;
while flag == 0
    if nargin < 4
        P = sumsub(V,B);
        if isempty(P)
            P = zeros(nx1,1);
        end
        Q = preimg(A,P);
        if isempty(Q)
            Q = zeros(nx1,1);
        end
        V = intsub(K,Q);
    else
        P = sumsub(V,B,tol);
        if isempty(P)
            P = zeros(nx1,1);
        end
        Q = preimg(A,P,tol);
        if isempty(Q)
            Q = zeros(nx1,1);
        end
        V = intsub(K,Q,tol);
    end
    m = size(V,2);
    if m == n
        flag = 1;
    else
        n = m;
    end  
end
W = V;
end