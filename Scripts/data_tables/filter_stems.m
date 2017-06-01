function  stems_filter = filter_stems( stems, tail_pos );
% stems_filter = filter_stems( stems, tail_pos );
%
% Filter stems to not include any with base pairs that extended into
%  flanking tails, as specified by tail_pos
%
% (C) R. Das, 2011, 2017

stems_filter = {};

for i = 1:length( stems )
  ok = 0;
  stem = stems{i};
  stem_new = [];
  
  for j = 1:size( stem, 1 )
    if ( isempty( find(stem(j,1) == tail_pos) )  & ...
	 isempty( find(stem(j,2) == tail_pos) )  )
      ok = 1;
      stem_new = [stem_new; stem(j,:) ];
    end
  end

  if ( ok ) 
    stems_filter = [ stems_filter, stem_new ]; 
  end;

end
