function  [ stem_correspondence, num_extra_bps_added_to_stem, num_missing_bps_in_stem ] = find_shared_stems( native_stems, stems );
%  [ stem_correspondence, num_extra_bps_added_to_stem, num_missing_bps_in_stem ] = find_shared_stems( native_stems, stems );
%
%  INPUTS:
%  native_stems = cell of Nnative arrays defining each native stem, each [Nbp x 2]
%  stems        = cell of Nmodel arrays defining each modeled stem, each [Nbp x 2]
%
%  OUTPUTS:
%  stem_correspondence = for each of Nnative native_stems, index of match in stems
%  num_extra_bps_added_to_stem = for each native_stem, num extra bps in
%                                    modeled match  (-1 means no match)
%  num_missing_bps_in_stem     = for each native_stem, num missing bps in
%                                    modeled match  (-1 means no match)
%

if nargin == 0;  help( mfilename ); return; end;

stem_correspondence          = zeros(  1,length( native_stems ) );
num_extra_bps_added_to_stem  = -1 * ones(  1,length( native_stems ) );
num_missing_bps_in_stem      = -1 * ones(  1,length( native_stems ) );

for i = 1:length( native_stems )
  native_stem = native_stems{i};
    
  for j = 1:length( stems )
    stem = stems{j};
  
    num_match = 0;

    for p = 1:size( native_stem, 1)
      found_match = 0;
      
      for q = 1:size( stem, 1)
	if ( native_stem( p, 1) == stem( q, 1 ) & ...
	     native_stem( p, 2) == stem( q, 2 ) )
	  found_match = 1; break;
	end
      end

      num_match = num_match + found_match;
            
    end
  
    if num_match
      stem_correspondence(i) = j;
      num_extra_bps_added_to_stem(i) = length(stem) - num_match;
      num_missing_bps_in_stem(i)     = length(native_stem) - num_match;
      break;
    end
    
  end
end

