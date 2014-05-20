EXE = 'python ../Utils/viennafold.py';
seq_file = 'bistable.seq';
pfs_file = 'bpp.txt';
command = [EXE,' ',seq_file,' ',pfs_file];
system( command );
bpp = load( 'bpp.txt' );
image( bpp*128 );
colormap( 1 - gray(100));
bpp(3,16)
bpp(16,27)
bp_ratio0 = bpp(3,16)/bpp( 16,27);
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NRES = 30;

SHAPE_file = 'SHAPE1.txt';
fid = fopen( SHAPE_file, 'w' );
SHAPE = -999 * ones(1,NRES);
SHAPE(3) = 0.0;
%SHAPE(26) = 0.0;
for i = 1:length(SHAPE); fprintf(fid, '%d %8.3f\n', i, SHAPE(i)); end;
fclose(fid );

command = [EXE,' ',seq_file,' ',pfs_file,' --sh ', SHAPE_file];
system( command );
bpp = load( 'bpp.txt' );
image( bpp*128 );
colormap( 1 - gray(100));
bpp(3,16)
bpp(16,27)
bp_ratio1 = bpp(3,16)/bpp(16,27);
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SHAPE_file = 'SHAPE2.txt';
fid = fopen( SHAPE_file, 'w' );
SHAPE = -999 * ones(1,NRES);
%SHAPE(26) = 0.3606;
SHAPE(27) = 0.0;
for i = 1:length(SHAPE); fprintf(fid, '%d %8.3f\n', i, SHAPE(i)); end;
fclose(fid );

command = [EXE,' ',seq_file,' ',pfs_file,' --sh ', SHAPE_file];
system( command );
bpp = load( 'bpp.txt' );
image( bpp*128 );
colormap( 1 - gray(100));
bpp(3,16)
bpp(16,27)
bp_ratio2 = bpp(3,16)/bpp( 16,27);
bp_ratio3 = bp_ratio1*bp_ratio2
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%expected = exp( -2 * 0.8/0.616);
expected = exp( -2 * 0.8/0.616);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EX_file = 'EX.txt';
fid = fopen( EX_file, 'w' );
EX = zeros( NRES,NRES );
EX( 3,16 ) = -0.80;
EX( 16,3 ) = -0.80;
%EX( 2,17 ) = -0.80;
%EX( 17,2 ) = -0.80;
for i = 1:NRES
  for j = 1:NRES
    fprintf(fid, '%8.3f,', EX(i,j)); 
  end;
  fprintf(fid, '\n');
end
fclose(fid );

command = [EXE,' ',seq_file,' ',pfs_file,' -x ', EX_file];
system( command );
bpp = load( 'bpp.txt' );
image( bpp*128 );
colormap( 1 - gray(100));
bpp(3,16)
bpp(16,27)
bp_ratio3 = bpp(3,16)/bpp( 16,27);

fprintf('\n\n');
fprintf( 'These should be equal:  favor_1/equal: %8.3f;   equal/favor_2  %8.3f\n',...
	 bp_ratio0/bp_ratio1, bp_ratio2/bp_ratio0 );
fprintf( 'And those should equal penalty for 2 * 0.8 kcal ==> %8.3f\n', expected );
fprintf( 'And here is the ratio for -ex pair bonus ==> %8.3f\n', bp_ratio0/bp_ratio3 );

