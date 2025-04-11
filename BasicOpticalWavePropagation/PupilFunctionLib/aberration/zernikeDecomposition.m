function [a,phs_recon] = zernikeDecomposition(phs_raw,nDecom,nflag)
%ZERNIKEDECOMPOSITION decomposes an phase matrix into different zernike
%modes
%
% Inputs-------------------------------------------------------------------
% phs_raw	-hybrid phase
% nDecom    -decompose 'phs_raw' into the first 'nDecom' zernike 
%            polynomials (Noll's index)
% nflag     -normalization flag
%
% Outputs------------------------------------------------------------------
% a         -amplitude of each zernike polynomial
% phs_recon -reconstructed phase
% phs_resid	-residual of between 'phs_raw' and 'phs_recon'
%
% Note---------------------------------------------------------------------
% This function can not be used for decomposite wrapped phase
%
% -------------------------------------------------------------------------
% Created by Xiang Hao,
% Modified by Xin Liu
% liuxin2018@zju.edu.cn
% Mar.19, 2021

if nargin==2
    nflag = [];
end

if size(phs_raw,1) ~= size(phs_raw,2)
    error('Input phase matrix mast be square!');
end

pv = max(phs_raw(:))-min(phs_raw(:));
if pv<=2*pi
    warning('The phase seems to have been wrapped, this function is not suitable!');
end

% pupil grid
[xx,yy] = meshgrid(linspace(-1,1,length(phs_raw)));
[theta,r] = cart2pol(xx,yy);
idx = (r<=1);

n = zeros(1,nDecom);
m = zeros(1,nDecom);
for ii = 1:nDecom
    [n(ii),m(ii)] = Noll2RA(ii);
end
zernModes = zernfun(n,m,r(idx),theta(idx),nflag);

a = zernModes\phs_raw(idx);  % least square solver of Z*A = PHI

phs_recon = zeros(size(phs_raw));
phs_recon(idx) = zernModes*a;  % recompute phase
end