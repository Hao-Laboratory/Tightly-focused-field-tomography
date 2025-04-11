function [amp,phs,varargout] = singleobjectivepr(Obj,PSF,Scope,pupilRes,Option)
%SINGLEOBJECTIVEPR Summary of this function goes here
%   Detailed explanation goes here
if ~isfield(Option,'UseGpu')
    Option.UseGpu = 0;
    warning(['UseGpu is set as default: ', num2str(Option.UseGpu)]);
end

if ~isfield(Option,'Precision')
    Option.Precision = 'double';
    warning(['Percision is set as default: ', Option.Precision]);
end

[amp,phs] = PR_MM(Obj,PSF,Scope,pupilRes,Option);
varargout{1} = [];
end