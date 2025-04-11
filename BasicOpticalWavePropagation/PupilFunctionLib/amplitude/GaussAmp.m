function amp = GaussAmp(pupilRes,w0)
% Gaussian amplitude distribution on the pupil plane.
%
% Input--------------------------------------------------------------------
% pupilRes: resolution of pupil function
% w0: relative width of beam waist
%
% -------------------------------------------------------------------------
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% May.15, 2020

[xp,yp] = meshgrid(linspace(-1,1,pupilRes));
[~, rho] = cart2pol(xp,yp);

amp = exp(-(rho/w0).^2);

amp(rho>1) = 0;
end