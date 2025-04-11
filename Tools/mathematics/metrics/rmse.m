function RMSE = rmse(M_sim,M_exp,varargin)
%RMSE calculates the Root-Mean-Square error between simulation and
%experimental data
%
% Input--------------------------------------------------------------------
% M_sim: simulation data (real or complex)
% M_exp: experiment data (real or complex)
%
% -------------------------------------------------------------------------
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Nov.19, 2020

switch nargin
    case 2
        normFlag = 1;
    case 3
        normFlag = varargin{1};
end

size1 = size(M_sim);
size2 = size(M_exp);

if ~isequal(size1,size2)
    error('Inputs should have the same dimensions!');
end

switch normFlag
    case 1
        RMSE = sqrt(sum(abs(M_sim-M_exp).^2,'all')./sum(abs(M_exp).^2,'all'));
    case 0
        RMSE = sqrt(sum(abs(M_sim-M_exp).^2,'all')./numel(M_exp));
    otherwise
        error('No such an option, please try again!');
end
end