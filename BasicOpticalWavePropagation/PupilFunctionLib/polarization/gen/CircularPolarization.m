function plr = CircularPolarization(pupilRes,handOption)
%CIRCULARPOLARIZATION generates circular state of polarization
%
% Author: Xin Liu
% email: liuxin2018@zju.edu.cn
% Mar.24, 2020

switch handOption
    case 'l'
        plr(:,:,1) = ones(pupilRes)./sqrt(2);
        plr(:,:,2) = ones(pupilRes)./sqrt(2).*1i;
        plr(:,:,3) = zeros(pupilRes);
        
    case 'r'
        plr(:,:,1) = ones(pupilRes)./sqrt(2);
        plr(:,:,2) = -ones(pupilRes)./sqrt(2).*1i;
        plr(:,:,3) = zeros(pupilRes);
    otherwise
        error('No such an option!');
end

x = linspace(-1,1,pupilRes);
[xx, yy] = meshgrid(x);
[~, rho] = cart2pol(xx,yy);

rho = repmat(rho,1,1,3);
plr(rho>1) = 0;
