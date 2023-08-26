function [R,V,F] = lctrb(K,A,B,tol,k)
% LCTRB Orthonormal basis computation for the largest controllability
% subspace of the pair (A,B) in the image of K. 
% [R,V,F] = lctrb(K,A,B) returns an orthonormal basis for the largest 
% controllability subspace of (A,B) in imK, an orthonormal basis for the
% largest (A,B)-controlled invariant subspace in  imK, and an F s.t. imV is
% (A+BF)-invariant. 
% 
% [R,V,F] = lctrb(K,A,B,tol) returns an orthonormal basis for the largest
% controllability subspace of (A,B) in imK, an orthonormal basis for the
% largest (A,B)-controlled invariant subspace in imK, and an F s.t. imV is 
% (A+BF)-invariant, where tol is a user-defined threshold to have a certain
% distance from the machine zero. This function uses isalci(K,A,B,tol) to 
% compute V and fci(V,A,B,tol) to compute F.
%
% [R,V,F] = lctrb(K,A,B,tol,k) returns an orthonormal basis for the largest
% controllability subspace of (A,B) in imK, an orthonormal basis for the 
% largest (A,B)-controlled invariant subspace in imK, and an F s.t. imV is 
% (A+BF)-invariant, where tol is a user-defined threshold to have a certain
% distance from the machine zero, and k is either 1 or 2. This function 
% uses isalci(K,A,B,tol) to compute V and fci(V,A,B,tol,k) to compute F.

% ========================================================================
% Selahattin Burak Sarsilmaz, February 2022 (Matlab R2021a), made use of 
% 1) Function BctrbL, Behcet Acikmese, September 2006. 
% Source: 
% i) Theorem 4.17 of Control theory for linear systems, 
% Trentelman, Stoorvogel, and Hautus. 
% ========================================================================

if isempty(A) || isempty(B) || isempty(K)
    error('The entries of lctrb must be nonempty matrices.');
end
[nx1,nx2] = size(A);
nx3 = size(B,1);
nv1 = size(K,1);
if nx1 ~= nx2 || nx2 ~= nx3 || nv1 ~= nx2
    error('Dimensions do NOT match.')
end

%1) Find the largest controlled invariant subspace for (A,B) in imK
if nargin > 3
    V = isalci(K,A,B,tol);
else
    V = isalci(K,A,B);
end
if isempty(V)
    Vne = zeros(nx1,1);
else
    Vne = V;
end

%2) Find a state feedback F s.t. imV is (A+BF)-invariant (F in F_(imV))
if nargin == 3
    F = fci(Vne,A,B);
end
if nargin == 4
    F = fci(Vne,A,B,tol);
end
if nargin == 5
    F = fci(Vne,A,B,tol,k);
end

%3) Find the intersection of imB and imV
if nargin > 3
    Bs = intsub(Vne,B,tol);
else
    Bs = intsub(Vne,B);
end

%4) Find the largest controllability subspace for (A,B) in imK
AF = A + B*F;
W = ctrb(AF,Bs);
if nargin > 3
    R = img(W,tol);
else
    R = img(W);
end
if isempty(R)
    Rne = zeros(nx1,1);
else
    Rne = R;
end

%5) Intersect R and V, and numerical warnings
if nargin > 3
    VintR = intsub(Vne,Rne,tol);
else
    VintR = intsub(Vne,Rne);
end
nv = size(V,2);
nr = size(R,2);
nvr = size(VintR,2);
if  nr ~= nvr
    error('Numerical issue: Use or update tol of lctrb.')
end
if  nv == nr && nvr == nv
    R = V;
end
end