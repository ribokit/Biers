[bpp_without, bpp_with] = function fold_with_ex_bonuses(sequence, EX)
fprintf('Folding without bonuses\n');
EXE = 'viennafold.py';
NRES = length(sequence)
s = [';\n',sequence,'\n1'];
fid = fopen('seqfile.seq', 'w');
fprintf(fid, s);
fclose(fid);
seq_file = 'seqfile.seq';
pfs_file = 'bpp.txt';
command = [EXE,' ',seq_file,' ',pfs_file];
system( command );
bpp_without = load( 'bpp.txt' );
figure(1)
title('Base pair probability matrix without bonuses')
image( bpp*128 );
colormap( 1 - gray(100));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Folding with bonuses\n');
fid = fopen( 'EX.txt', 'w' );
for i = 1:NRES
  for j = 1:NRES
    fprintf(fid, '%8.3f,', EX(i,j)); 
  end;
  fprintf(fid, '\n');
end
fclose(fid );

command = [EXE,' ',seq_file,' ',pfs_file,' -x EX.txt'];
system( command );
bpp_with = load( 'bpp.txt' );
figure(2)
title('Base pair probabilities with bonuses')
image( bpp*128 );
colormap( 1 - gray(100));
end
