function [structure_mohca, bpp_2D_mohca, d_out, structure_na] = mohca_RNAstructure(name, rdatname, structure, outprefix, pseudoknots, nboot, zsthreshold, docombine, dobroaden, binarize)

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
    d_diag_combined = combine_diagonal(d_final, 0, 0.025, 2);
    d_final = d_diag_combined;
else
    d_diag_combined = [];
end
if( dobroaden )
    d_broadened = broaden(d_final, 1);
    d_final = d_broadened;
else
    d_broadened = [];
end
if( binarize )
    d_final(d_final > 0) = 1;
end

seqstr = {sequence, structure};

mohcaplot_biers(d_final, '', '', [outprefix, name, ': Secondary Structure Bonuses from MOHCA'], '', [outprefix, name, '_bonuses.eps'], seqstr, 1 );

if( zsthreshold < 0 )
    d_final = [];
else
    d_final = -d_final(1:length(sequence),1:length(sequence));
end

[structure_mohca, bpp_2D_mohca] = rna_structure( sequence, [], 0, [], d_final, nboot, pseudoknots );
structure_na = rna_structure( sequence, [], 0, [], [], 0, pseudoknots );
fprintf(1, '%s\n', structure_mohca);
fprintf(1, '%s\n', structure);
if isempty(structure); 
    str1 = structure_mohca; 
    str2 = structure_na;
else
    str1 = structure;
    str2 = structure;
end;

output_varna_html( [outprefix, name, '_MOHCA.html'], sequence, structure_mohca, str1, structure_mohca, rdat.offset, [], [], [], bpp_2D_mohca );
output_varna_html( [outprefix, name, '_NA.html'], sequence, structure_na, str2, structure_na, rdat.offset, [], [], [], [] );

mohcaplot_biers(-d_final, '', '', [outprefix, name, ': NoData Overlay'], '', [outprefix, name, '_NA_overlay.eps'], {sequence, structure_na}, 2 );
mohcaplot_biers(-d_final, '', '', [outprefix, name, ': MOHCA Overlay'], '', [outprefix, name, '_MOHCA_overlay.eps'], {sequence, structure_mohca}, 2 );

d_out = { d_final, d_filtered, d_diag_combined, d_broadened };
end