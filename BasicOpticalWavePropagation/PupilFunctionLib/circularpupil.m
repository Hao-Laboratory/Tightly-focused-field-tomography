function pupil_arg = circularpupil(pupil_arg)
%CIRCULARPUPIL add circular pupil to input matrix
%
% Xin Liu
% liuxin2018@zju.edu.cn
% Nov.15, 2021

M = length(pupil_arg);
x = linspace(-1,1,M);
y = x;
[xx,yy] = meshgrid(x,y);
rho = xx.^2 + yy.^2;

rho = repmat(rho,1,1,size(pupil_arg,3));
pupil_arg(rho>1) = 0;
end