function [Vmin, V, F] = ssbci(K,A,B,E,tol,k)
% SSBCI Orthonormal basis computation for the smallest self bounded
% (A,B)-controlled invariant subspace contained in the image of K and 
% containing the image of E.
% [Vmin,V,F] = ssbci(K,A,B,E) returns an orthonormal basis for the smallest
% self bounded (A,B)-controlled invariant subspace contained in imK and 
% containing imE, an orthonormal basis for the largest (A,B)-controlled 
% invariant subspace in imK, and an F s.t. imV is (A+BF)-invariant. 
%
% [Vmin,V,F] = ssbci(K,A,B,E,tol) returns an orthonormal basis for the 
% smallest self bounded (A,B)-controlled invariant subspace contained in 
% imK and containing imE, an orthonormal basis for the largest 
% (A,B)-controlled invariant subspace in imK, and an F s.t. imV is 
% (A+BF)-invariant, where tol is a user-defined threshold to have a certain 
% distance from the machine zero. This function uses isalci(K,A,B,tol) to 
% compute V and fci(V,A,B,tol) to compute F.
%
% [Vmin,V,F] = ssbci(K,A,B,E,tol,k) returns an orthonormal basis for the 
% smallest self bounded (A,B)-controlled invariant subspace contained in 
% imK and containing imE, an orthonormal basis for the largest 
% (A,B)-controlled invariant subspace in imK, and an F s.t. imV is 
% (A+BF)-invariant, where tol is a user-defined threshold to have a certain
% distance from the machine zero, and k is either 1 or 2. This function 
% uses isalci(K,A,B,tol) to compute V and fci(V,A,B,tol,k) to compute F.

% ========================================================================
% Selahattin Burak Sarsilmaz, February 2022 (Matlab R2021a).
% Source: 
% i) Multivariable regulation in geometric terms: Old and new results,
% Marro.
% ========================================================================

if isempty(A) || isempty(B) || isempty(K) || isempty(E)
    error('The entries of ssbci must be nonempty matrices.');
end
[nx1,nx2] = size(A);
nx3 = size(B,1);
nv1 = size(K,1);
ne1 = size(E,1);
if nx1 ~= nx2 || nx2 ~= nx3 || nv1 ~= nx2 || ne1 ~= nx2
    error('Dimensions do NOT match.')
end

%1) Find the largest (A,B)-invariant subspace for (A,B) in imK
if nargin > 4
    V = isalci(K,A,B,tol);
else
    V = isalci(K,A,B);
end
if isempty(V)
    Vne = zeros(nx1,1);
else
    Vne = V;
end

% 2) Check whether the set of (A,B)-invariant subspaces contained in imK
% and containing imE is nonempty
if nargin > 4
    E = img(E,tol);
else
    E = img(E); 
end
if isempty(E)
    Ene = zeros(nx1,1);
else
    Ene = E;
end
if nargin > 4
    VintE = intsub(Ene,Vne,tol);
else
    VintE = intsub(Ene,Vne);
end
nv = size(V,2);
ne = size(E,2);
nve = size(VintE,2);
if  ne~= nve 
    error('EITHER there exists no (A,B)-controlled invariant subspace contained in and containing the given subspaces OR there is a numerical issue: Use or update tol of ssbci.')
end
if ne == nv && nve == ne
    V = E; 
end
if isempty(V)
    Vne = zeros(nx1,1);
else
    Vne = V;
end

% 3) Find a state feedback F s.t. imV is (A+BF)-invariant (F in F_(imV))
if nargin == 4
    F = fci(Vne,A,B);
end
if nargin == 5
    F = fci(Vne,A,B,tol);
end
if nargin == 6
    F = fci(Vne,A,B,tol,k);
end

%4) Find the intersection of imV and (imB + imE)
if nargin > 4
    BsE = sumsub(B,E,tol);
    if isempty(BsE)
        BsEne = zeros(nx3,1);
    else
        BsEne = BsE;
    end
    Bs = intsub(Vne,BsEne,tol);
else
    BsE = sumsub(B,E);
    if isempty(BsE)
        BsEne = zeros(nx3,1);
    else
        BsEne = BsE;
    end
    Bs = intsub(Vne,BsEne);
end

% 5) Find the smallest self-bounded (A,B)-invariant subspace contained in
% imK and containing imE
AF = A + B*F;
W = ctrb(AF,Bs);
if nargin > 4
    Vmin = img(W,tol);
else
    Vmin = img(W);
end
if isempty(Vmin)
    Vminne = zeros(nx1,1);
else
    Vminne = Vmin;
end

% 6) Intersect Vmin and V, and numerical warnings
if nargin > 4
    VintVmin = intsub(Vne,Vminne,tol);
else
    VintVmin = intsub(Vne,Vminne);
end
nvmin = size(Vmin,2);
nvvmin = size(VintVmin,2);
if  nvmin ~= nvvmin
    error('Numerical issue: Use or update tol in ssbci.')
end
if  nv == nvmin && nvvmin == nv
    Vmin = V;
end
end