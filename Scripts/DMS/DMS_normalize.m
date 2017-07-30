function [d_norm, scalefactor, cap_value ] = DMS_normalize( d, sequence );
% [d_norm, scalefactor, cap_value ] = DMS_normalize( d, sequence );
%
%  'box-plot' normalize for DMS data:
%  -- set values at any positions that are not A or C to NaN so they're not considered for
%     normalization.
%  -- remove outliers, i.e., any values above 1.5 * interquartile range 
%  -- Find maximum value after filtering. If its smaller than the 95th percentile value, then 
%      take that instead.
%  -- scalefactor is mean of top 10th percentile of values, but removing values above that filter.
%
% INPUTS:
%   d        = reactivity
%   sequence = sequence; must be same length as d.
%
% OUTPUTS:
%   d_norm   = box-plot normalized reactivity
% scalefactor= how much d was scaled to get d_norm
% cap_value  = cutoff for outliers.
%
% Note that this is basically a wrapper around SHAPE_normalize. 
% 
% (C) Das lab, Stanford University, 2017.

if nargin < 2;  help( mfilename ); return; end;
if length( sequence ) ~= size( d, 1 ); error( 'Size of sequence must match d' ); end;

DMS_pos = [strfind(sequence,'A'),strfind( sequence,'C')]; %,strfind(r.sequence,'a'),strfind(r.sequence,'c')];
non_DMS_pos = setdiff( [1:size( d, 1 )], DMS_pos );
d( non_DMS_pos, : ) = NaN;
[ d_norm, scalefactor, cap_value ] = SHAPE_normalize( d );
