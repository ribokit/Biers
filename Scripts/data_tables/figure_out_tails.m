function  tail_pos = figure_out_tails( structure );

tail_pos = [];
i = 1;
while structure(i) == '.'; tail_pos = [ tail_pos, i ]; i = i+1; end
i = length(structure);
while structure(i) == '.'; tail_pos = [ tail_pos, i ]; i = i-1; end
tail_pos = sort( tail_pos );

