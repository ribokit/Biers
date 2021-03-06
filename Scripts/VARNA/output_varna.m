function [special_base_pairs, special_colors] = output_varna(varna_file, sequence, display_structure, native_structure, modeled_structure, offset, ROI, crystpos, shape, bpp, color_mode, flat)
% [special_base_pairs, special_colors] = output_varna( ...
%                               varna_file, sequence, display_structure, native_structure, modeled_structure, 
%                               offset, ROI, crystpos, shape, bpp, color_mode); 
%
% Generate HTML file with VARNA applet for RNA secstr visualziation.
% **Make sure to set your path in get_varna_exe.m
%
% [Input]
% varna_file            Required            File name (must end in: .html, .eps, .svg, .png)                                                          .png
% sequence              Required            RNA sequence for display
% display_structure     Required            RNA secstr for display
% native_structure      Optional            Reference/native secstr used for secstr comparison.
% modeled_structure     Optional            Predicted/new secstr used for secstr comparison.
% offset                Optional            Sequence numbering offset, default 0.
% ROI                   Optional            Region of interest (don't show
%                                                base pairs outside this range)
% crystpos              Optional            [?]
% shape                 Optional            SHAPE reactivity profile for nucleotide coloring, default none.
% bpp                   Optional            Base-pairing probability matrix for helix-wise confidence score.
% color_mode            Optional            Secstr comparison line color set, default 1.
% flat                  Optional            make exterior loop flat, default: 1.
%
% [Output]
% special_base_pairs    Index of base-pairs that differ between native_structure and modeled_structure
% special_colors        Secstr comparison color set
%
% (C) T47, Stanford University 2013-2015
% (C) Rhiju Das, Stanford University 2017

if nargin == 0;  help( mfilename ); return; end;

if ~exist('native_structure','var') | isempty(native_structure); native_structure = display_structure; end;
if ~exist('modeled_structure','var') | isempty(modeled_structure); modeled_structure = display_structure; end;
if ~exist('bpp', 'var') | isempty(bpp); bpp = []; end;
if ~exist('shape', 'var') | isempty(shape); shape = []; end;

if ~exist('offset','var') | isempty(offset); offset = 0; end;
if ~exist('ROI','var') | isempty(ROI); ROI = [1:length(sequence)] + offset; end;
if ~exist('crystpos','var') | isempty(crystpos); crystpos = [1:length(sequence)] + offset; end;
if ~exist('color_mode', 'var') | isempty(color_mode); color_mode = 1; end;
if ~exist('flat', 'var') | isempty(flat); flat = 1; end;

ROI = sort(ROI);
seqpos = [1:length(sequence)] + offset;

% select subset of residues based on ROI
gp = find(seqpos >= min(ROI)  &  seqpos <= max(ROI));
ROI_offset = min(gp) - 1;
crystpos_subset = sort(crystpos) - offset - ROI_offset;

% Need to shift around numbers if we have only focused on certain residues (ROI).
sequence_subset = sequence(gp);
display_structure_subset = get_subset_structure(display_structure, ROI - offset);
%native_structure_subset  = native_structure;
native_structure_subset  = get_subset_structure(native_structure, ROI - offset);
modeled_structure_subset = get_subset_structure(modeled_structure, ROI - offset);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define specially colored base pairs.
special_base_pairs = {}; special_colors = {};

[added_base_pairs, missing_base_pairs, added_base_pairs_EDGE, missing_base_pairs_EDGE] = compare_structures(modeled_structure_subset, native_structure_subset);
  
not_in_cryst_base_pairs = [];
[added_base_pairs, not_in_cryst_base_pairs] = filter_crystpos(added_base_pairs, not_in_cryst_base_pairs, crystpos_subset);
[added_base_pairs_EDGE, not_in_cryst_base_pairs] = filter_crystpos(added_base_pairs_EDGE, not_in_cryst_base_pairs, crystpos_subset);

special_base_pairs{1} = missing_base_pairs;
special_base_pairs{2} = missing_base_pairs_EDGE;
special_base_pairs{3} = added_base_pairs;
special_base_pairs{4} = added_base_pairs_EDGE;

if color_mode == 1;
  special_colors{1} = [1, 0, 1];
  special_colors{2} = [1, 0, 1];
  special_colors{3} = [0, 1, 1];
  special_colors{4} = [0, 1, 1];
elseif color_mode == 2;
  special_colors{1} = [1, 0.5, 0];
  special_colors{2} = [1, 0.5, 0];
  special_colors{3} = [0.5, 0.5, 0.5];
  special_colors{4} = [0.5, 0.5, 0.5];
else
  special_colors{1} = [0, 0.5, 0];
  special_colors{2} = [0, 0.5, 0];
  special_colors{3} = [0, 0.5, 1];
  special_colors{4} = [0, 0.5, 1];
end;


if ~isempty(not_in_cryst_base_pairs);
  special_base_pairs{5} = not_in_cryst_base_pairs;
  special_colors{5} = [0.7, 0.7, 0.7];
end;

bpp_values = [];
bpp_anchor_bases = [];
if ~isempty(bpp);
  bpp_subset = bpp(gp, gp);
  stems = parse_stems(display_structure_subset);
  [bpp_values, bpp_anchor_bases] = get_bpp_values(stems, bpp_subset); 
end;


% Do it!
fprintf( ['Creating file for VARNA visualization: ', varna_file, '\n' ] );

shape_subset = [];
if ~isempty(shape); 
    for i = length(shape)+1: length(sequence); shape(i) = NaN; end; % fill with NaNs
    shape_subset = shape(gp); 
end;

ADJUST_SHAPE = 1.0;
varna_fig(varna_file, sequence_subset, display_structure_subset, ...
    ADJUST_SHAPE * shape_subset, 2 , offset + ROI_offset, ...
    special_base_pairs, special_colors, bpp_values, bpp_anchor_bases, flat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function structure_shift = get_subset_structure(structure, ROI)

bps = convert_structure_to_bps(structure);
bps_shift = [];

for i = 1:size(bps, 1);
  bp_shift1 = find(bps(i, 1) == ROI);
  bp_shift2 = find(bps(i, 2) == ROI);
  if (~isempty(bp_shift1) && ~isempty(bp_shift2));
    bps_shift = [bps_shift; bp_shift1, bp_shift2];
  end;
end;

structure_shift = convert_bps_to_structure(bps_shift, length(ROI));
structure_shift = structure;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [added_base_pairs_filter, not_in_cryst_base_pairs ] = filter_crystpos( added_base_pairs, not_in_cryst_base_pairs, crystpos_subset )

added_base_pairs_filter = [];
for i = 1:size(added_base_pairs, 1);
  if ( isempty( find( added_base_pairs(i,1) == crystpos_subset ) ) || ...
       isempty( find( added_base_pairs(i,2) == crystpos_subset ) ) )
    not_in_cryst_base_pairs = [not_in_cryst_base_pairs; added_base_pairs(i, :)];
  else
    added_base_pairs_filter = [added_base_pairs_filter; added_base_pairs(i, :)];
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [ bpp_values, anchor_bases ] = get_bpp_values(stems, bpp_subset)
bpp_values = [];
anchor_bases = [];
for i = 1:length(stems);
  stem = stems{i};
  middle_index = round(length( stem ) / 2);
  anchor_bases(i) = stem(middle_index, 1);
  
  bpp_value = 0.0;
  for j = 1: size(stem, 1);
    bpp_value = max(bpp_value,  bpp_subset(stem(j, 1), stem(j, 2)));
  end;
  bpp_values(i) = bpp_value;
end;

