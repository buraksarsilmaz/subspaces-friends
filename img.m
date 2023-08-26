function W = img(X,tol,p)
% IMG  Orthonormal basis computation for the image of a matrix.
% W = img(X) returns an orthonormal basis for the image of X based on
% the built-in function 'orth'. 
%
% W = img(X,tol) returns an orthonormal basis for the image of X based on 
% the built-in function 'svd', where tol is a user-defined threshold to
% have a certain distance from the machine zero.
% 
% W = img(X,tol,1) returns an orthonormal basis for the image of X based on
% the built-in function 'qr', where tol is a user-defined threshold to
% have a certain distance from the machine zero.
% 
% W = img(X,tol,2) returns an ordered orthonormal basis for the image of X 
% based on the built-in function 'qr', where tol is a user-defined 
% threshold to have a certain distance from the machine zero.

% ========================================================================
% Selahattin Burak Sarsilmaz, February 2022 (Matlab R2021a), made use of 
% 1) Function Bim, Behcet Acikmese,  June 2007;
% 2) Function ima, Giuseppe Basile & Giovanni Marro , May 1997.
% ========================================================================

[m,n] = size(X);
%% Use orth, which relies on svd
if nargin == 1
     W = orth(X);
     if size(W,2) == m
         W = eye(m);
     end
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
    if r > 0
        W = U(:,1:r);
        if  r == m
            W = eye(m);
        end
    else
        W = zeros(m,0);
    end
end
%% Use qr with tolerance
if nargin == 3 
    if p == 1
    [Q,R,P]=qr(X); 
    if m == 1 || n == 1
       d = abs(R(1,1));      
    else
       d = abs(diag(R));  
    end
    nz = find(d>tol);
    r = length(nz);    % Rank of X
        if r > 0
        W = Q(:,1:r);
            if  r == m
             W = eye(m);
            end
        else
        W = zeros(m,0);
        end
    end
    if p == 2
    i = 1;
    X1 = X;
        while i == 1
            na = size(X1,2);
            v = 1:na;
            [Q,R]=qr(X1);  
            if m == 1 || na == 1
                d = abs(R(1,1)); 
            else
                d = abs(diag(R)); 
            end
            fz = find(d<=tol,1);
            if ~isempty(fz)
                v = find(v ~= fz);
            end
            if isempty(fz) || isempty(v)
               i = 0;
            else 
               X1 = X1(:,v);
            end
        end
        nz = find(d>tol);
        r = length(nz);    % Rank of X
        if r > 0
            W = Q (:,1:r);
        else 
            W = zeros(m,0);
        end 
    end  
end
end