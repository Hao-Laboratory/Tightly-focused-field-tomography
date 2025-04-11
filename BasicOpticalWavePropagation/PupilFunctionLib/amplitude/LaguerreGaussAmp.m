function amp = LaguerreGaussAmp(pupilRes,w0,l)
% LAGUERREGAUSSAMP generates the Laguerre-Gauss [LG(0,l) modes]amplitude distribution on 
% the pupil plane.
%
% Input--------------------------------------------------------------------
% pupilRes: resolution of pupil function
% w0: beam waist
% l: the topological charge
%
% -------------------------------------------------------------------------
% Author: Liu, Xin
% Email:liuxin2018@zju.edu.cn
% Mar.14, 2020

x = linspace(-1,1,pupilRes);
[xx, yy] = meshgrid(x);
[~, rho] = cart2pol(xx,yy);

tr = rho/w0;
amp = (sqrt(2).*tr).^abs(l).*exp(-tr.^2);

amp(rho>=1) = 0;
end