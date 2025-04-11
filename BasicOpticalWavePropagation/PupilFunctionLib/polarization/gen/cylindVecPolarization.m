function plr = cylindVecPolarization(pupilRes,theta,n)
%VECPOLARIZATION generates cylinderical vector beams, theta is the angle 
% between the polarization vector and the radial vector.
%
% Input--------------------------------------------------------------------
% n > 1 for high-order vector beam
% theta = 0 or pi: radial polarization
% theta = pi/2 or -pi/2: azimuthal polarization
% theta = pi/4*0.935, with uniform amplitude can form flat-top PSF
% -------------------------------------------------------------------------
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Jun.09, 2020
% -------------------------------------------------------------------------

x = linspace(-1,1,pupilRes);
[xx, yy] = meshgrid(x);
[phi, rho] = cart2pol(xx,yy);

plr(:,:,1) = cos(n*phi+theta);
plr(:,:,2) = sin(n*phi+theta);
plr(:,:,3) = zeros(pupilRes);

rho = repmat(rho,1,1,3);
plr(rho>1) = 0;
end