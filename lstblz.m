function Vs = lstblz(K,A,B,tol,k)
% LSTBLZ Orthonormal basis computation for the largest stabilizability
% subspace of the pair (A,B) in the image of K. 
% Vs = lstblz(K,A,B) returns an orthonormal basis for the largest 
% stabilizability subspace of the pair (A,B) contained in imK. 
%
% Vs = lstblz(K,A,B,tol) returns an orthonormal basis for the largest 
% stabilizability subspace of the pair (A,B) contained in imK, where 
% tol is a user-defined threshold to have a certain distance from the
% machine zero. This function uses lctrb(K,A,B,tol).  
% 
% Vs = lstblz(K,A,B,tol,k) returns an orthonormal basis for the largest 
% stabilizability subspace of the pair (A,B) contained in imK, where 
% tol is a user-defined threshold to have a certain distance from the
% machine zero, and k is either 1 or 2. This function uses 
% lctrb(K,A,B,tol,k).  

% ========================================================================
% Selahattin Burak Sarsilmaz, August 2023 (Matlab R2023a), made use of 
% 1) Function BstabL, Behcet Acikmese, June 2007. 
% Source: 
% i) Corollary 4.27 of Control theory for linear systems, 
% Trentelman, Stoorvogel, and Hautus. 
% ========================================================================

if isempty(A) || isempty(B) || isempty(K)
    error('The entries of lstblz must be nonempty matrices.');
end
[nx1,nx2] = size(A);
nx3 = size(B,1);
nv1 = size(K,1);
if nx1 ~= nx2 || nx2 ~= nx3 || nv1 ~= nx2
    error('Dimensions do NOT match.')
end

% Find the largest controllability subspace imR for (A,B) in imK, the 
% largest (A,B)-controlled invariant subspace imV in imK, and a 
% state feedback F s.t. imV is (A+BF)-invariant
if nargin == 3
    [R, V, F]= lctrb(K,A,B); 
end
if nargin == 4
    [R, V, F]= lctrb(K,A,B,tol); 
end
if nargin == 5
    [R, V, F]= lctrb(K,A,B,tol,k); 
end 

if isempty(V)
    Vne = zeros(nx1,1);
else
    Vne = V;
end

AF = A+B*F; 

% Find the stable subpsace of AF and the largest stabilizability subspace
% contained in imK. 
if nargin > 3
    X = stabsub(AF,tol);
    if isempty(X)
        Xne = zeros(nx1,1);
    else
        Xne = X;
    end
    Int = intsub(Xne,Vne,tol);
    Vs = sumsub(Int,R,tol);
else 
    X = stabsub(AF);
    if isempty(X)
        Xne = zeros(nx1,1);
    else
        Xne = X;
    end
    Int = intsub(Xne,Vne);
    Vs = sumsub(Int,R);
end

if isempty(Vs)
    Vsne = zeros(nx1,1);
else
    Vsne = Vs;
end

% Intersect Vs and V, and numerical warnings
if nargin > 3
    VintVs = intsub(Vne,Vsne,tol);
else
    VintVs = intsub(Vne,Vsne);
end
nv = size(V,2);
ns = size(Vs,2);
nvs = size(VintVs,2);
if  ns ~= nvs
    error('Numerical issue: Use or update tol of lstblz.')
end
if  nv == ns && nvs == nv
    Vs = V;
end

end