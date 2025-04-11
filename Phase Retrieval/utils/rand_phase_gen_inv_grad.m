function rand_phs = rand_phase_gen_inv_grad(Obj,Beam,roi_size_x,roi_size_y,rand_Level)
%RAND_PHASE_GEN_INV_GRAD generate random phase with controlled gradient
%using inverse gradient methods
%   Detailed explanation goes here
%
% Xin Liu
% liuxin2018@zju.edu.cn
% Nov.22, 2023

% rand_Level = 5;

kc = 2*pi*Obj.NA/Beam.wavelength;
dk = 2*kc/(Beam.PupilRes-1);

x = linspace(-1,1,rand_Level); y = x; [xx,yy] = meshgrid(x,y);
xq = linspace(-1,1,Beam.PupilRes); yq = xq; [xxq,yyq] = meshgrid(xq,yq);
[~,rho] = cart2pol(xxq,yyq);
mask = rho<=1;

phs_x_ini = rand(rand_Level,rand_Level);
phs_y_ini = rand(rand_Level,rand_Level);

phs_grad_x_tmp = interp2(xx,yy,phs_x_ini,xxq,yyq,'spline');
phs_grad_y_tmp = interp2(xx,yy,phs_y_ini,xxq,yyq,'spline');

phs_grad_x = zeros(Beam.PupilRes);
phs_grad_y = zeros(Beam.PupilRes);
phs_grad_x(mask) = rescale(phs_grad_x_tmp(mask),-roi_size_x,roi_size_x);
phs_grad_y(mask) = rescale(phs_grad_y_tmp(mask),-roi_size_y,roi_size_y);

rand_phs = intgrad2(phs_grad_x,phs_grad_y,dk,dk);

rand_phs(~mask) = 0;
rand_phs(mask) = rand_phs(mask) - min(rand_phs(mask));

% rand_phs = wrapTo2Pi(rand_phs);

% verification
% figure;
% [FX,FY] = gradient(rand_phs,dk,dk);
% FX(rho>0.97) = 0;
% FY(rho>0.97) = 0;
% 
% subplot(2,2,1);imagesc(phs_grad_x);axis image xy off;colorbar;title('dphi_x GT');
% subplot(2,2,2);imagesc(phs_grad_y);axis image xy off;colorbar;title('dphi_y GT');
% subplot(2,2,3);imagesc(FX);axis image xy off;colorbar;title('dphi_x gen');
% subplot(2,2,4);imagesc(FY);axis image xy off;colorbar;title('dphi_y gen');
end
