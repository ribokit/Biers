function rZ = get_Zscore_rdat( r, Z )
% rZ = get_Zscore_rdat( r, Z )
%
% Repackages RDAT after Z-score transformation, useful for 
%  preparing input to map2dplot
%
% INPUTS:
% r  = RDAT object for a mutate-and-map (or M2-seq) data set
% Z  = Z-score matrix derived from, e.g., output_Zscore_from_rdat
% 
% OUTPUTS:
% rZ = output RDAT object 
%
% (C) R. Das, Stanford University,  2017
%

rZ = r;
rZ.reactivity = []; % zero out first row.
rZ.reactivity(:,[2:size(r.reactivity,2)]) = -Z;
rZ.annotations = [ r.annotations, 'processing:Zscore' ];
