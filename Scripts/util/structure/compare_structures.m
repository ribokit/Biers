function  [added_base_pairs, missing_base_pairs, added_base_pairs_EDGE, missing_base_pairs_EDGE ] = compare_structures( structure, native_structure );

bps        = convert_structure_to_bps( structure );
native_bps = convert_structure_to_bps( native_structure );

[ missing_base_pairs, missing_base_pairs_EDGE]   = get_missing_base_pairs( native_bps, bps );
[   added_base_pairs,  added_base_pairs_EDGE ]  = get_missing_base_pairs( bps, native_bps );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ missing_base_pairs, missing_base_pairs_EDGE ] = get_missing_base_pairs( bps, native_bps );

missing_base_pairs = [];
missing_base_pairs_EDGE = [];

for i = 1:size( bps, 1 );

  if ~check_bp( native_bps, bps(i,:) )
    missing_base_pairs = [ missing_base_pairs; bps(i,:) ];
    
    % adjacent to and stacked onto native base pair
    if check_bp( native_bps, bps(i,:) + [-1 +1] ) 
      missing_base_pairs_EDGE = [missing_base_pairs_EDGE; bps(i,:) ];
    end
    if check_bp( native_bps, bps(i,:) + [+1 -1] ) 
      missing_base_pairs_EDGE = [missing_base_pairs_EDGE; bps(i,:) ];
    end
    if check_bp( bps, bps(i,:) + [-1 +1] ) & ...
	  check_bp( native_bps, bps(i,:) + [-2 +2] )
      missing_base_pairs_EDGE = [missing_base_pairs_EDGE; bps(i,:) ];
    end
    if check_bp( bps, bps(i,:) + [+1 -1] ) & ...
	  check_bp( native_bps, bps(i,:) + [+2 -2] )
      missing_base_pairs_EDGE = [missing_base_pairs_EDGE; bps(i,:) ];
    end

  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = check_bp( native_bps, bp )

ok = 1;
if length( native_bps ) == 0; ok = 0; return; end;

idx = find( bp(1) == native_bps(:,1) );

if isempty( idx ) | ( bp(2) ~= native_bps(idx,2) )
  ok = 0;
end    
