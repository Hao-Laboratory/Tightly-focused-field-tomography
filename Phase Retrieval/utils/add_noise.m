function noisy_img = add_noise(clear_img,std)
%ADD_NOISE add shot noise and dark/readout noise to optical images
%   
% clear_img: unit with the number of photons
% std: level (standard deviation) of Gaussian noise, unit with the number of photons
% 
% Xin Liu
% liuxin2018@zju.edu.cn
% Dec. 18, 2023

shot_noise_img = poissrnd(clear_img);
dark_readout_noise = std*randn(size(clear_img));
noisy_img =  shot_noise_img + dark_readout_noise;
noisy_img(noisy_img<0) = 0;
end

