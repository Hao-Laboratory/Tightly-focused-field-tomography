function Beam = pupilreini(Beam,Option)
% PUPILREINI reinitialize pupil function with random distribution
%
% Xin Liu
% liuxin.optics@gmail.com

if contains(Option.PrMode,'amp')
    Beam.amp = ones(Beam.PupilRes);
    
    if strcmp(Option.Precision,'single')
        Beam.amp = single(Beam.amp);
    end
    if Option.UseGpu == 1
        Beam.amp = gpuArray(Beam.amp);
    end
end

if contains(Option.PrMode,'phs')
    Beam.phs = zeros(Beam.PupilRes);
    Beam.abr = zeros(Beam.PupilRes);
    if strcmp(Option.Precision,'single')
        Beam.phs = single(Beam.phs);
        Beam.abr = single(Beam.abr);
    end
    if Option.UseGpu == 1
        Beam.phs = gpuArray(Beam.phs);
        Beam.abr = gpuArray(Beam.abr);
    end
else
    Beam.phs = Beam.phs + Beam.abr;
end

if contains(Option.PrMode,'plr')
    Beam.plr = randomPolarization(Beam.PupilRes,4);
    if strcmp(Option.Precision,'single')
        Beam.plr = single(Beam.plr);
    end
    if Option.UseGpu == 1
        Beam.plr = gpuArray(Beam.plr);
    end
end

end