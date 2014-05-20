function run_test_case(name, rdatname, structure, outprefix, zsthreshold, docombine, dobroaden, binarize)

rdat = read_rdat_file(rdatname);
sequence = rdat.sequence
xpos = strfind(sequence, 'X');
sequence =  sequence(1:xpos(2)-1);
sequence = strrep(sequence, 'X', '');
d = mohca_symmetric(rdat.reactivity);
d_filtered = diagonal_filter(d, zsthreshold);
d_final = d_filtered;
if( docombine )
    d_final = combine_diagonal(d_final, 0, 0.025, 2);
end
if( dobroaden )
    d_final = broaden(d_final, 1);
end
if( binarize )
    d_final(d_final > 0) = 1;
end

plot_and_save_sec_str_MOHCA(d_final, sequence, structure, [outprefix, name, '_bonuses.eps'])

if( zsthreshold < 0 )
    d_final = [];
else
    d_final = -d_final(1:length(sequence),1:length(sequence));
end

[structure_mohca, bpp_2D] = rna_structure( sequence, [], 0, [], d_final, 0, 1 );
fprintf(1, '%s\n', structure_mohca);
fprintf(1, '%s\n', structure);
output_varna_html( [outprefix, name, '.html'], sequence, structure_mohca, structure, structure_mohca, rdat.offset, [], [], [], bpp_2D );


end