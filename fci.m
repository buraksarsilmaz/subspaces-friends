function F = fci(V,A,B,tol,k)
% FCI  A state feedback F such that image of V is (A+BF)-invariant. 
% F = fci(V,A,B) returns an F such that imV is (A+BF)-invariant iff
% imV is (A,B)-controlled invariant, where two systems of linear equations 
% are solved consecutively via the built-in function 'pinv'.
%
% F = fci(V,A,B,tol) returns an F such that imV is (A+BF)-invariant iff
% imV is (A,B)-controlled invariant, where two systems of linear equations 
% are solved via the built-in function 'pinv' and tol is a user-defined 
% threshold to have a certain distance from the machine zero.
% 
% F = fci(V,A,B,tol,k) returns an F such that imV is (A+BF)-invariant iff
% imV is (A,B)-controlled invariant, where one system of linear equations 
% is solved, tol is a user-defined threshold to have a certain distance 
% from the machine zero, and k is either 1 or 2. When k = 1, the built-in 
% function 'mldivide' is used. When k = 2, the built-in function 'pinv' is 
% used. 

% ========================================================================
% Selahattin Burak Sarsilmaz, February 2022 (Matlab R2021a), made use of 
% 1) Function BciF, Behcet Acikmese,  June 2007. 
% Sources: 
% i) Theorem 4.2 of Control theory for linear systems, 
% Trentelman, Stoorvogel, and Hautus;
% ii) Remark 2.2 of Revisiting disturbance decoupling with optimization 
% perspective, Sarsilmaz, Li, and Acikmese.
% ========================================================================

if isempty(A) || isempty(B) || isempty(V)
    error('The entries of fci must be nonempty matrices.');
end
[nx1,nx2] = size(A);
nx3 = size(B,1);
nv1 = size(V,1);
if nx1 ~= nx2 || nx2 ~= nx3 || nv1 ~= nx2
    error('Dimensions do NOT match.')
end
 
if nargin == 3
    V = img(V);
    nv = size(V,2);
    rhs = A*V;
    AV = img(rhs);
    n = size(AV,2);
    if isempty(AV)
        AV = zeros(nx1,1);
    end
    VpB = sumsub(V,B);
    if isempty(VpB)
        VpB = zeros(nx1,1);
    end
    W = intsub(AV,VpB);
    m = size(W,2);
else
    V = img(V,tol);
    nv = size(V,2);
    rhs = A*V;
    AV = img(rhs,tol);
    n = size(AV,2);
    if isempty(AV)
        AV = zeros(nx1,1);
    end
    VpB = sumsub(V,B,tol);
    if isempty(VpB)
        VpB = zeros(nx1,1);
    end
    W = intsub(AV,VpB,tol);
    m = size(W,2);
end

q = size(B,2);
if n ~= m
    error('EITHER given subspace is not (A,B)-controlled invariant OR there is a numerical issue: Use or update tol of fci.');
else
    if nargin <= 4
        lhs = [V B];
        %lhs x = rhs  has a solution iff imV is controlled invariant.
        %[V B] * [X' U']' = A*V, where U = -F*V.
        sol = pinv(lhs)*rhs;
        U = sol(nv+1:end,:);
        Ft = pinv(V')*(-U');
        F = Ft';
    else
        % V*X - B*F*V = A*V
        vrhs = rhs(:);
        vlhs = [kron(eye(nv),V) (-1)*kron(V',B)];
        if k == 1
            sol = vlhs\vrhs; % underdetermined system 
        end
        if k == 2
            % least-squares solution
            sol = pinv(vlhs)*vrhs;
        end
        vF = sol(nv^2+1:end);
        F = reshape(vF,q,nv1);
    end
end
end