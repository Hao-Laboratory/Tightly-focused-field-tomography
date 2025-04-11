function mat_out = angle2pi(mat_in)
%ANGLE2PI return the angle of input in 2pi radian
%
% LIU Xin
% liuxin24@hku.hk
% May 26, 2024

mat_out = wrapTo2Pi(angle(mat_in));
end

