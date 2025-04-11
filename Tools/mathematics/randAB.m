function M = randAB(a,b,sz)
%RANDAB generate random matrix between in the range of a~b, where a<b
%   Detailed explanation goes here

% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Feb.23, 2021

M = a + (b-a)*rand(sz);
end