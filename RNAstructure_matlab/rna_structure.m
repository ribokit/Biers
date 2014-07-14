function [structure, bpp, SHAPE_out ] = rna_structure( sequence, area_shape, offset, seqpos, EX, NUM_BOOTSTRAP, cmd_pk, area_dms, temperature, shape_intercept, shape_slope, maxdist )
% [structure, bpp, SHAPE_out ] = rna_structure( sequence, area_shape, offset, seqpos, EX, NUM_BOOTSTRAP, cmd...
%                                               temperature, shape_intercept, shape_slope, maxdist);
%
% make sure to set your path in get_exe_dir.m
%

if ~exist( 'area_shape','var' ); area_shape = []; end;
if ~exist( 'area_dms','var' ); area_dms = []; end;
if ~exist( 'offset','var' ); offset = 0; end;
if ~exist( 'seqpos','var' ); seqpos = 0; end;
if ~exist( 'EX','var' ); EX = []; end;
if ~exist( 'NUM_BOOTSTRAP','var' ); NUM_BOOTSTRAP = 0; end;
if ~exist( 'cmd_pk','var') || isempty(cmd_pk); 
    cmd_pk = 0; 
elseif cmd_pk == 1;
    fprintf('WARNING: ShapeKnots run is much slower than Fold!\n');
end;
if ~exist( 'maxdist','var' ); maxdist = 0; end;
seq_file = 'tmp.seq';
fid = fopen( seq_file, 'w' );
fprintf( fid, ';\n' );
fprintf( fid, 'default\n' );
fprintf( fid, '%s1\n',sequence );
fclose( fid );

area_shape( find( isnan( area_shape ) ) ) = -999;
area_dms( find( isnan( area_dms ) ) ) = -999;

bps = [];
nres = length( sequence );

EX_file = '';
if ~isempty( EX )
  EX_file = 'tmp_EX.txt';
  save 'tmp_EX.txt' -ascii EX;
end


SHAPE_out = [];

if ~isempty( area_shape )
  area_shape_norm = area_shape; %SHAPE_normalize( area_shape );
  clf
  plot( seqpos, area_shape_norm ); 
  
  SHAPE_file = ['tmp_SHAPE.txt'];
  fid = fopen( SHAPE_file, 'w' );

  SHAPE_out = nan * ones( 1, nres );
  for i = 1:nres
    gp =  find( i + offset == seqpos );
    if ~isempty( gp )
      fprintf( fid, '%d  %8.3f\n', i, area_shape_norm( gp ) );
      SHAPE_out( i ) = area_shape_norm( gp );
    end
  end
  fclose( fid );
else
  SHAPE_file = '';
end

if ~isempty( area_dms )
  area_dms_norm = area_dms; %SHAPE_normalize( area_shape );
  figure();clf;
  plot( seqpos, area_dms_norm ); 
  
  DMS_file = ['tmp_DMS.txt'];
  fid = fopen( DMS_file, 'w' );

  DMS_out = nan * ones( 1, nres );
  for i = 1:nres
    gp =  find( i + offset == seqpos );
    if ~isempty( gp )
      fprintf( fid, '%d  %8.3f\n', i, area_dms_norm( gp ) );
      DMS_out( i ) = area_dms_norm( gp );
    end
  end
  fclose( fid );
else
  DMS_file = '';
end

if ~exist( 'temperature','var' ); temperature = 24; end;

if ~exist( 'shape_intercept','var' ); shape_intercepts = []; 
else  shape_intercepts = [ shape_intercept ]; end;

if ~exist( 'shape_slope','var' ); shape_slopes = []; 
else  shape_slopes = [ shape_slope ]; end;

experimental_offset = 0.0;
zscore_scalings = [];
USE_VIENNA = 0;

tic;
[bpp, structure, ct_file, command ] = run_rna_structure_with_bootstrap( NUM_BOOTSTRAP, seq_file, bps, offset, nres, EX_file, SHAPE_file, ...
    temperature, experimental_offset, zscore_scalings, shape_intercepts, shape_slopes, USE_VIENNA, maxdist, cmd_pk, DMS_file);
toc;

% beep notice when finished
gong;


