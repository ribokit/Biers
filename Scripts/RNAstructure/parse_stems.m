function stems = parse_stems( structure )
% stems = parse_stems( structure )
%
%  INPUT
%   structure = structure in dot parens notation
% 
%  OUTPUT
%   stems = cell of stems, each as an array of base pairs.
%
bps   = convert_structure_to_bps( structure );
stems = parse_stems_from_bps( bps );

