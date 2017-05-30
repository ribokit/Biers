function bps = convert_structure_to_bps( structure );
%  bps = convert_structure_to_bps( structure );
%
%  INPUTS
%  structure = structure in dot-parens notation, using characters .()[]{}
%  
%  OUTPUT
%  bps  = matrix of base pairs.
%
% (C) Rhiju Das, Stanford University, 2011-2012, 2017

bps = [];
bps = get_bps( structure, bps, '(', ')' );
bps = get_bps( structure, bps, '[', ']' );
bps = get_bps( structure, bps, '{', '}' );

function bps = get_bps( structure, bps, left_delim, right_delim );

LEFT_BRACKET = [];
for i = 1:length(structure )
  switch structure(i)
   case left_delim
    LEFT_BRACKET = [LEFT_BRACKET, i];
   case right_delim
    bps = [bps; LEFT_BRACKET(end), i];
    LEFT_BRACKET = LEFT_BRACKET( 1 : end-1 );
  end
end

