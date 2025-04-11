function [amp,phs,plr] = pupilfuntion_decode(Epx,Epy)
%PUPILFUNTION_DECODE decode amplitude, phase, and polarization from two
%orthognal polarization state of pupil function
%   
% LIU Xin
% liuxin2018@zju.edu.cn
% Jun.18, 2023

% pupil amplitude
pupil_amp = sqrt(abs(Epx).^2+abs(Epy).^2);

% pupil phase
pupil_phs = angle(Epx+Epy);

% pupil polarization
absXp = abs(Epx)./pupil_amp;
absYp = abs(Epy)./pupil_amp;
angleXp = angle(Epx)-pupil_phs;
angleYp = angle(Epy)-pupil_phs;
pupil_px = absXp.*exp(1i*angleXp);
pupil_py = absYp.*exp(1i*angleYp);

% amp
amp = pupil_amp;
amp(~isfinite(amp)) = 0;

% phase
phs = pupil_phs;
phs(~isfinite(phs)) = 0;

% plr
pupil_px(~isfinite(pupil_px)) = 0;
pupil_py(~isfinite(pupil_py)) = 0;
plr(:,:,1) = pupil_px;
plr(:,:,2) = pupil_py;
end

