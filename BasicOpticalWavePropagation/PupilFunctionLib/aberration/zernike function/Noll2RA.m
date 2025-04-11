% NATURALINDEX transfer the Noll's indices (j) of Zernike polynomials to 
% the natural arrange of the indices n (radial index) and m (azimuthal).
%
% format: [n,m] = Noll2RA(j)
%
% See also NollIndex
% coded by HAO, Xiang
% xiang.hao@yale.edu
% Jan 20, 2017

function [n,m] = Noll2RA(j)

%% calculate n
counter = 1;
n = 1;

while j > counter
    n = n + 1;
    counter = counter + n;
end
n = n - 1;

%% calculate m
rawSequence = j - n * ( n + 1) / 2;

if mod(n,2) == 0 % n is even
    m = ceil((rawSequence - 1) / 2) * 2;
else             % n is odd
    m = floor((rawSequence - 1) / 2) * 2 + 1;
end

% odd js are assigned to m<0
% if mod(j,2) == 1
if j - floor(j/2)*2 == 1
    m = -m;
end

end