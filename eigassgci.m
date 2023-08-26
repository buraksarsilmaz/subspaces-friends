function [fxdeigXmodS, fxdeigVmodR, pr, psv, F] = eigassgci(V,A,B,tol,k)
% EIGASSGCI Computes the fixed eigenvalues of the given (A,B)-controlled 
% invariant subspace imV and assigns all the assignable eigenvalues to 
% desired locations under the constraint that imV is (A+BF)-invariant. 
% [fxdeigXmodS, fxdeigVmodR] =  eigassgci(V,A,B,tol) produces the fixed
% eigenvalues of the linear maps (A+BF)|imX/imS and (A+BF)|imV/imR for all
% F in F(imV), respectively, where imX is the set of n-dimensional real 
% vectors, imS = imV + <A,imB>, imR is the the largest controllability 
% subspace for (A,B) in imV, and tol is a user-defined threshold to have a 
% certain distance from the machine zero. This function uses lctrb(K,A,B,tol).
%
% [fxdeigXmodS, fxdeigVmodR, pr, psv, F] =  eigassgci(V,A,B,tol) also 
% produces state feedback F in F(imV) such that all the assignable 
% eigenvalues are assigned to desired locations pr and psv, where pr and 
% psv are the user-defined eigenvalues of the linear maps (A+BF)|imR and 
% (A+BF)|imS/imV, respectively.
%
% [fxdeigXmodS, fxdeigVmodR] =  eigassgci(V,A,B,tol,k) produces the fixed
% eigenvalues of the linear maps (A+BF)|imX/imS and (A+BF)|imV/imR for all
% F in F(imV), respectively, where imX is the set of n-dimensional real 
% vectors, imS = imV + <A,imB>, imR is the the largest controllability 
% subspace for (A,B) in imV, tol is a user-defined threshold to have a 
% certain distance from the machine zero, and k is either 1 or 2. This 
% function uses lctrb(K,A,B,tol,k).
%
% [fxdeigXmodS, fxdeigVmodR, pr, psv, F] =  eigassgci(V,A,B,tol,k) also 
% produces state feedback F in F(imV) such that all the assignable 
% eigenvalues are assigned to desired locations pr and psv, where pr and 
% psv are the user-defined eigenvalues of the linear maps (A+BF)|imR and
% (A+BF)|imS/imV, respectively. 

% ========================================================================
% Selahattin Burak Sarsilmaz, August 2023 (Matlab R2023a) 
% Source: 
% i) The constructive proof of Theorem 4.18 in Control theory for 
% linear systems, Trentelman, Stoorvogel, and Hautus. 
% ========================================================================

if isempty(A) || isempty(B) || isempty(V)
    error('The entries of eigassgci must be nonempty matrices.');
end
[nx1,nx2] = size(A);
nx3 = size(B,1);
nv1 = size(V,1);
if nx1 ~= nx2 || nx2 ~= nx3 || nv1 ~= nx2
    error('Dimensions do NOT match.')
end

% 1) Find a state feedback F0 s.t. imV is (A+BF0)-invariant (F0 in F_(imV))
if nargin > 4
    F0 = fci(V,A,B,tol,k);
else
    F0 = fci(V,A,B,tol);
end

% 2) Find the largest controllability subspace for (A,B) in imV
if nargin > 4 
     R = lctrb(V,A,B,tol,k);
else
     R = lctrb(V,A,B,tol);
end

% 3) Compute the strongly invariant subspace S = imV + <A,imB>
V = img(V,tol);
S = sumsub(V,ctrb(A,B),tol);

% 4) Find a basis of nx2-dimensional Euclidean space adapted to the chain:
% imR subset imV subset imS subset nx2-dimensional Euclidean space
r = size(R,2);
v = size(V,2);
s = size(S,2); 
P_pre = [R V S eye(nx2)];
P = img(P_pre,tol,2); 

% 5) Algebraically equivalent system
Atil = P\(A*P); 
Btil = P\B;

% 6) Fixed eigenvalues of the map restricted to nx2-dimensional Euclidean 
% space mod imS
fxdeigXmodS = []; 
if s < nx2
    Atil44 = Atil(s+1:nx2,s+1:nx2);
    fxdeigXmodS = esort(eig(Atil44));
end

% 7) Fixed eigenvalues of the map restricted to imV mod imR for all F in 
% F(imV)
fxdeigVmodR = [];
AF0 = A + B*F0;
AF0til = P\(AF0*P); 
if r < v
    AF0til22 = AF0til(r+1:v,r+1:v);
    fxdeigVmodR = esort(eig(AF0til22));
    nfu = find(real(fxdeigVmodR) >= 0,1);
    if isempty(nfu)
        disp('V is internally stabilizable.')
    else 
        disp('V is NOT internally stabilizable.')
    end
else 
    disp('V is internally stabilizable.')
end
F0til = F0*P;

% 8) Assign the assignable eigenvalues to desired locations under the 
% constraint that F is in F(imV)
if nargout > 2
    pr =[];
    psv =[];
    if r > 0
        AF0til11 = AF0til(1:r,1:r);
        fprintf('specify %i',r)
        fprintf(' eigenvalue(s) for (A+BF)|R\n')
        pr = input('in a scalar or vector: '); 
        L = preimg(B,V,tol);
        k = size(L,2); 
        BtilL = Btil*L; 
        Btil1L = BtilL(1:r,k);
        F11til = -place(AF0til11,Btil1L,pr);
        F1til = [F11til zeros(k,nx2-r)]; 
        F2til = F0til + L*F1til; 
    else
        F2til = F0til;
    end
    if s-v > 0
        fprintf('specify %i',s-v)
        fprintf(' eigenvalue(s) for (A+BF)| S/V\n')
        psv = input('in a scalar or vector: ');
        AF2til = Atil + Btil*F2til; 
        AF2til33 = AF2til(v+1:s,v+1:s);
        Btil3 = Btil(v+1:s,:);
        Ftil33 = -place(AF2til33,Btil3,psv);
        m = size(B,2);
        F3til = [zeros(m,v) Ftil33 zeros(m,nx2-s)];
        Ftil = F2til +F3til; 
    else
        Ftil = F2til;
    end
    F = Ftil/P;
end

end