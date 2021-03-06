function [sens, ppv] = calc_accuracy_bp(structures, native_str)

if ~exist('native_str','var') || isempty(native_str); return; end;

sens = zeros(1,length(structures));
ppv = zeros(1,length(structures));

for i = 1:length(structures);
    [sens(i), ppv(i)] = calc_sens_ppv_bp(structures{i}, native_str);
end;
