function conventional_stem_colors = get_conventional_stem_colors( native_stems );
% conventional_stem_colors = get_conventional_stem_colors( native_stems );
%
%  assigns stems into 'domains' based on whether they are connected by
%  short internal loops (2-way junctions). Then assigns some pretty colors.
%
% (C) R. Das, 2017.

% Colors for each stem
% cred, ... are attempt to reproduce "clarence" colors from M2-seq P4-P6
% depictions.
cred = [1 0.2 0.5]; ctan = [1,0.5,0.3]; cblue = [0.1 0.5 0.8]; cgold = [0.9, 0.5 0.3]; 
cgreen = [0.2 0.7 0.2]; cpurple = [1 0.2 1]; cpink = [1 0.8 0.8]; ccyan = [0.1 0.5 1];
conventional_colors = [cred;ctan;cblue;cgold;cgreen;cpurple; cpink; ccyan];
% fill out with MATLAB colors.
conventional_colors = [conventional_colors; jet(length(native_stems) - length(conventional_colors))];

conventional_stem_colors = conventional_colors( define_domains( native_stems ), : );
