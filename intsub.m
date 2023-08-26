function W = intsub(X,Y,tol)
% INTSUB Orthonormal basis computation for the intersection of two subspaces
% W = intsub(X,Y) returns an orthonormal basis for the intersection of two
% subspaces: imX and imY, via default 'kern' and 'sumsub'.
%
% W = intsub(X,Y,tol) returns an orthonormal basis for the intersection of
% two subspaces: imX and imY via 'kern' and 'sumsub', where tol is a 
% user-defined threshold to have a certain distance from the machine zero.

% ========================================================================
% Selahattin Burak Sarsilmaz, February 2022 (Matlab R2021a), made use of 
% 1) Function BintS, Behcet Acikmese,  September 2006;
% ========================================================================

if isempty(X) || isempty(Y) 
    error('The entries of intsub must be nonempty matrices.');
end
m1 = size(X,1);
m2 = size(Y,1);
if m1 ~= m2
    error('Subspaces do NOT lie in the same vector space.');
end
%% Use default kern and sumsub
if nargin == 2
% imX intersect imY = (kerX'+kerY')perb = (imV1+imV2)perb = (imVs)perb 
% = kerVs', where imVs = imV1 + imV2
V1 = kern(X');
V2 = kern(Y');
Vs = sumsub(V1,V2);
W = kern(Vs');
end
%% Use svd-based kern and sumsub with tol
if nargin == 3
V1 = kern(X',tol);
V2 = kern(Y',tol);
Vs = sumsub(V1,V2,tol);
W = kern(Vs',tol);
end
end