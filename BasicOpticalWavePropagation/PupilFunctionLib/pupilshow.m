function pupilshow(pupil_arg,circularFlag)
%PUPILSHOW show amplitude and phase of pupil inside a circular aperture
%(set the values out of pupil to be transparent)
%
% Xin Liu
% liuxin2018@zju.edu.cn
% Nov.16, 2020

if nargin == 1
    circularFlag = 1;
end
M = length(pupil_arg);
x = linspace(-1,1,M);
y = x;
[xx,yy] = meshgrid(x,y);
r = xx.^2 + yy.^2;

alpha = ones(M);
if circularFlag == 1
    alpha(r>1) = 0;
end
imagesc(pupil_arg,'AlphaData',alpha);
axis image xy off;
end