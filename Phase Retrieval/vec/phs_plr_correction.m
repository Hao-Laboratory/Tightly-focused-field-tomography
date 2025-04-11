function [phsC,plrC] = phs_plr_correction(plr_ini,plr_res,phs_res)
%PHS_PLR_CORRECTION correct the phase and polarization of phase retrieval
%result
%   
% INPUT********************************************************************
% plr_ini: expected polarization of pupil function
% plr_res: restored polarization of pupil function
% phs_res: restored phase of pupil function
%
% OUTPUT*******************************************************************
% phsC: corrected phase compared with plr_ini
% plrC: corrected polarization which is similar to plr_ini
%
% *************************************************************************
% LIU Xin
% liuxin2018@zju.edu.cn
% Apr.26, 2021

plrC = zeros(size(plr_res));

psiX = angle(plr_res(:,:,1).*conj(plr_ini(:,:,1)));
psiY = angle(plr_res(:,:,2).*conj(plr_ini(:,:,2)));
psi = angle(exp(1i*psiX)+exp(1i*psiY));

% corrected phase
phsC = wrapTo2Pi(phs_res+psi);

% corrected polarization
plrC(:,:,1) = plr_res(:,:,1).*exp(-1i*psi);
plrC(:,:,2) = plr_res(:,:,2).*exp(-1i*psi);
end

