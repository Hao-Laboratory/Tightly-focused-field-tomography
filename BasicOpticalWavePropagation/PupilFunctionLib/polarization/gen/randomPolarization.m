function plr = randomPolarization(pupilRes,RandLevel)
%RANDOMPOLARPOLARIZATION generates random polarization state
%   
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Jul.12, 2021
xy = linspace(-1,1,pupilRes);
[xx, yy] = meshgrid(xy);
[~, rho] = cart2pol(xx,yy);

% aviod too small amplitude
minVal = 0.2;
maxVal = 0.8;
% minVal = 0;
% maxVal = 1;
Seed = imresize(randAB(minVal,maxVal,RandLevel),[pupilRes pupilRes],'bicubic');
Seed = clamp(Seed,minVal,maxVal);
Ax0 = sqrt(Seed);
Ay0 = sqrt(1-Ax0.^2);

deltaX = 2*pi.*imresize(rand(RandLevel),[pupilRes pupilRes],'bicubic');
deltaY = 2*pi.*imresize(rand(RandLevel),[pupilRes pupilRes],'bicubic');

plr(:,:,1) = Ax0.*exp(1i*deltaX);
plr(:,:,2) = Ay0.*exp(1i*deltaY);
plr(:,:,3) = zeros(pupilRes);

rho = repmat(rho,1,1,3);
plr(rho>1) = 0;
end

