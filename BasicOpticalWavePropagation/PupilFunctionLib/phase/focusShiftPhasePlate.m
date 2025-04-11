function phs = focusShiftPhasePlate(pupilRes,wavelength,NA,n,dx,dy,dz)
%FOCUSSHIFTPHASEPLATE shifts the focus in the image space accurately.
% 
% this function is based on the paper "Dual Coaxial Longitudinal 
% Polarization Vortex Structures"
%
% pupilDiaPixNum: the number of pixels in the pupil.
% wavelength: the specific wavelength for defocus.
% NA: the numerical aperture of the objective lens.
% n: the RI in image space.
%
% dz: the magnitude of defocus. 
% dz>0 means the focus shifts along positive direction of z-axis. 
% dz<0 means the focus shifts along negative direction of z-axis. 
%
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Mar.13, 2020

x = linspace(-1,1,pupilRes);
[xx, yy] = meshgrid(x);
[phi, rho] = cart2pol(xx,yy);

sinTheta = NA/n.*rho;
sinTheta(rho>1) = 0;
theta = asin(sinTheta);

k = 2*pi*n/wavelength;
kx = -k.*sin(theta).*cos(phi);
ky = -k.*sin(theta).*sin(phi);
kz = k.*cos(theta);

phs = (-kx.*dx-ky.*dy-kz.*dz);

phs(rho>1) = 0;
% phs = mod(phs,2*pi);
end
