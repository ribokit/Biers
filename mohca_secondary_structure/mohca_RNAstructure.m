function mohca_RNAstructure(name, rdatname, structure, outprefix, pseudoknots, nboot, zsthreshold, docombine, dobroaden, binarize)

if ~exist('zsthreshold', 'var')
    zsthreshold = 1.2;
end
if ~exist('outprefix', 'var')
    outprefix = '';
end
if ~exist('docombine', 'var')
    docombine = true;
end
if ~exist('dobroaden', 'var')
    dobroaden = true;
end
if ~exist('binarize', 'var')
    binarize = true;
end
if ~exist('nboot', 'var')
    nboot = 0;
end
if ~exist('pseudoknots', 'var')
    pseudoknots = 0;
end

rdat = read_rdat_file(rdatname);
sequence = rdat.sequence;
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

seqstr = {sequence, structure};

mohcaplot_biers(d_final, '', '', [outprefix, name, ': Secondary Structure Bonuses from MOHCA'], '', [outprefix, name, '_bonuses.eps'], seqstr );

if( zsthreshold < 0 )
    d_final = [];
else
    d_final = -d_final(1:length(sequence),1:length(sequence));
end

[structure_mohca, bpp_2D] = rna_structure( sequence, [], 0, [], d_final, nboot, pseudoknots );
fprintf(1, '%s\n', structure_mohca);
fprintf(1, '%s\n', structure);
output_varna_html( [outprefix, name, '.html'], sequence, structure_mohca, structure, structure_mohca, rdat.offset, [], [], [], bpp_2D );

end