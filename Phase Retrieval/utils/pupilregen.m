function beam = pupilregen(beam,Epx,Epy,sys)
%PUPILREGEN regenerate pupil function (amplitude, phase and polarization)
%from Epx and Epy
%
% Xin Liu
% liuxin2018@zju.edu.cn
% May.8, 2023

% pupil amplitude
pupil_amp = sqrt(abs(Epx).^2+abs(Epy).^2);

% pupil phase
if contains(sys.prMode,'plr')
    pupil_phs = angle(Epx+Epy);
%     pupil_phs = angle(Epx.^2+Epy.^2);
%     pupil_phs = angle(sqrt(Epx.^2+Epy.^2));
    %         pupil_phs = angle(Epx);
    %         pupil_phs = angle(mean(cat(3,Epx,Epy),3));
%     [Ep_p,Ep_s] = Exyz2Eps(Epx,Epy,0,0,sys.phi);
%     pupil_phs = angle(Ep_p+Ep_s);
    
else
    pupil_phs = angle(Epx.*conj(beam.plr(:,:,1)) +...
        Epy.*conj(beam.plr(:,:,2)));
    %         pupil_phs = angle(mean(cat(3,Epx.*conj(beam.plr(:,:,1)),...
    %             Epy.*conj(beam.plr(:,:,2))),3));
end

% pupil polarization
absXp = abs(Epx)./pupil_amp;
absYp = abs(Epy)./pupil_amp;
if contains(sys.prMode,'phs')
    angleXp = angle(Epx)-pupil_phs;
    angleYp = angle(Epy)-pupil_phs;
else
    angleXp = angle(Epx)-beam.phs;
    angleYp = angle(Epy)-beam.phs;
end
pupil_px = absXp.*exp(1i*angleXp);
pupil_py = absYp.*exp(1i*angleYp);

if contains(sys.prMode,'amp')
    pupil_amp(sys.rho>1) = 0;
    beam.amp = pupil_amp;
    beam.amp(~isfinite(beam.amp)) = 0;
end

if contains(sys.prMode,'phs')
    pupil_phs(sys.rho>1) = 0;
    beam.phs = pupil_phs;
    beam.phs(~isfinite(beam.phs)) = 0;
end

if contains(sys.prMode,'plr')
    pupil_px(~isfinite(pupil_px)) = 0;
    pupil_py(~isfinite(pupil_py)) = 0;
    beam.plr(:,:,1) = pupil_px;
    beam.plr(:,:,2) = pupil_py;
    beam.plr(:,:,3) = zeros(size(pupil_px));
end
end

