function [amp,phs] = simple_PR_MM(sys,PHI,PSF)
% simple_PR_MM is the simplified version of 'PR_MM.m'

% electric field distribution in focal space
Ef = PSF.amp.*exp(1i*PSF.phs);

EHold = pagemtimes(pagemtimes(sys.My',Ef),sys.Mx');

%% pupil function extraction 
Ep = EHold.*conj(PHI);
Ep = mean(Ep,3);  % this averages the gradient

% pupil function
Ep = Ep./sys.prefix;

% energy conservation
Ep = Ep.*sys.dx*sys.dy;

amp = abs(Ep);
phs = angle(Ep);
end