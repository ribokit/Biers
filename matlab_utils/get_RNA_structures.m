function [structures, energies] = get_RNA_structures(sequence, number, bonuses)
FOLD = '/home/tsuname/Documents/lab/rhiju/rdat/external/RNAstructure/exe/Fold';
CT2DOT = '/home/tsuname/Documents/lab/rhiju/rdat/external/RNAstructure/exe/ct2dot';
seq_file = '/tmp/rnastructure.seq';
ct_file = '/tmp/rnastructure.ct';
shape_file = '/tmp/shape.txt';
seqfid = fopen(seq_file, 'w');
fprintf(seqfid, ';\nblah\n%s1', sequence);
fclose(seqfid);
if isempty(bonuses)
    command = [FOLD,' ',seq_file,' ',ct_file, ' -sh ', shape_file];
else 
    shpfid = fopen(shape_file, 'w');
    for i=1:length(bonuses)
        if ~isnan(bonuses(i))
          fprintf(shpfid, '%d %8.3f\n', i, bonuses(i));
        end
    end
    fclose(shpfid);
    command = [FOLD,' ',seq_file,' ',ct_file];
end
structures = {};
energies = {};
system( command );
for j=1:number
    system( [CT2DOT, ' ', ct_file, ' ', num2str(j), ' /tmp/rnastructure.dot']);
    dotfid = fopen('/tmp/rnastructure.dot');
    tline = fgetl(dotfid);
    while ischar(tline)
        if tline(1) == ')' || tline(1) == '.' || tline(1) == '('
            structures{j} = tline;
        elseif tline(1) == '>'
	    energies{j} = tline;
	end
        tline = fgetl(dotfid);
    end
    fclose(dotfid);
end
end
