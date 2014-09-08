function  [ special_base_pairs, special_colors ] = output_varna_html( varna_file, sequence, display_structure, native_structure, modeled_structure, offset, mutpos, crystpos, shape, bpp )

if ~exist( 'mutpos','var' ) || isempty( mutpos);  mutpos = [1:length(sequence)]+offset; end;
if ~exist( 'crystpos','var' ) || isempty(crystpos);  crystpos = [1:length(sequence)]+offset; end;

mutpos = sort( mutpos );
seqpos = [1:length(sequence)] + offset;

% select subset of residues based on mutpos
gp = find(  seqpos >= min( mutpos )  &  seqpos <= max(mutpos ) );
mutpos_offset = min(gp) - 1;
crystpos_subset = sort(crystpos) - offset - mutpos_offset;

% Need to shift around numbers if we have only focused on certain residues (mutpos).
sequence_subset = sequence( gp );
display_structure_subset = get_subset_structure( display_structure, mutpos-offset );
%native_structure_subset  = native_structure;
native_structure_subset  = get_subset_structure( native_structure, mutpos-offset );
modeled_structure_subset = get_subset_structure( modeled_structure, mutpos-offset );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define specially colored base pairs.
special_base_pairs = {}; special_colors = {};

[added_base_pairs, missing_base_pairs, added_base_pairs_EDGE, missing_base_pairs_EDGE ] = ...
    compare_structures( modeled_structure_subset, native_structure_subset );
  
not_in_cryst_base_pairs = [];
[added_base_pairs, not_in_cryst_base_pairs ] = filter_crystpos( added_base_pairs, not_in_cryst_base_pairs, crystpos_subset );
[added_base_pairs_EDGE, not_in_cryst_base_pairs ] = filter_crystpos( added_base_pairs_EDGE, not_in_cryst_base_pairs, crystpos_subset );

special_base_pairs{1} = missing_base_pairs;
special_colors{1} = [1,0.5,0];
%special_colors{1} = [0,0.5,0];

special_base_pairs{2} = missing_base_pairs_EDGE;
special_colors{2} = [1.0,0.5,0];
%special_colors{2} = [0,0.5,0];

special_base_pairs{3} = added_base_pairs;
%special_colors{3} = [0.5,0.5,0.5];
special_colors{3} = [0,0.5,1];

special_base_pairs{4} = added_base_pairs_EDGE;
special_colors{4} = [0,0.5,1];

if ~isempty( not_in_cryst_base_pairs )
  special_base_pairs{5} = not_in_cryst_base_pairs;
  special_colors{5} = [0.7, 0.7, 0.7 ];
end

bpp_values = [];
bpp_anchor_bases = [];
if ~isempty( bpp )
  bpp_subset = bpp( gp, gp );
  stems = parse_stems( display_structure_subset );
  [ bpp_values, bpp_anchor_bases ] = get_bpp_values( stems, bpp_subset ); 
end


% Do it!
fprintf( ['Creating html file for VARNA visualization: ', varna_file, '\n' ] );

shape_subset = [];
if ~isempty( shape ); shape_subset = shape( gp ); end;

ADJUST_SHAPE = 1.0;
varna_fig( varna_file,sequence_subset,display_structure_subset,ADJUST_SHAPE*shape_subset,2, offset + mutpos_offset, special_base_pairs, special_colors, bpp_values, bpp_anchor_bases); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function structure_shift = get_subset_structure( structure, mutpos )

bps = convert_structure_to_bps( structure );
bps_shift = [];

for i = 1:size( bps,1 )
  bp_shift1 = find(bps(i,1) == mutpos);
  bp_shift2 = find(bps(i,2) == mutpos);
  if ( ~isempty( bp_shift1 ) && ~isempty( bp_shift2 ) ) 
    bps_shift = [ bps_shift; bp_shift1, bp_shift2 ];
  end
end

structure_shift = convert_bps_to_structure( bps_shift, length(mutpos) );
structure_shift = structure;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [added_base_pairs_filter, not_in_cryst_base_pairs ] = filter_crystpos( added_base_pairs, not_in_cryst_base_pairs, crystpos_subset )

added_base_pairs_filter = [];
for i = 1:size( added_base_pairs, 1 )
  if ( isempty( find( added_base_pairs(i,1) == crystpos_subset ) ) || ...
       isempty( find( added_base_pairs(i,2) == crystpos_subset ) ) )
    not_in_cryst_base_pairs = [ not_in_cryst_base_pairs; added_base_pairs(i,:) ];
  else
    added_base_pairs_filter = [ added_base_pairs_filter; added_base_pairs(i,:) ];
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [ bpp_values, anchor_bases ] = get_bpp_values( stems, bpp_subset )

for i = 1:length( stems );
  stem = stems{i};
  middle_index = round( length( stem ) / 2 );
  anchor_bases(i) = stem( middle_index, 1 );
  
  bpp_value = 0.0;
  for j = 1: size( stem, 1 )
    bpp_value = max( bpp_value,  bpp_subset( stem(j,1), stem(j,2) ) );
  end
  bpp_values(i) = bpp_value;
end

