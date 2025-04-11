function mask = genCircularMask(pupilRes)
%GENCIRCULARMASK generates a circular mask
%
% LIU Xin
% liuxin2018@zju.edu.cn
% Apr. 22, 2022

r0 = (pupilRes-1)/2;
x = linspace(-r0,r0,pupilRes);
y = x;
[xx,yy] = meshgrid(x,y);
r = hypot(xx,yy);

mask = true(pupilRes);
mask(r>r0) = false;
end

