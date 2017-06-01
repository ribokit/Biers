function [ crystpos, tail_crystpos ] = get_tail_crystpos( tag, structure )

outpath1 = 'OUTPUT1D/';
crystpos_file = [outpath1,tag,'_crystpos.txt'];
crystpos = textread( crystpos_file, '%d' );
offset_file = [outpath1,tag,'_offset.txt'];
offset = textread( offset_file, '%d' ); 
tail_crystpos = setdiff( [1:length(structure)], crystpos-offset );
