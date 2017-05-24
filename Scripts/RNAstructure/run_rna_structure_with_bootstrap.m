function [bpp, structure, ct_file, command ] = run_rna_structure_with_bootstrap( NUM_BOOTSTRAP, seq_file, bps, offset, nres, EX_file, SHAPE_file, ...
    temperature, experimental_offset, zscore_scaling, shape_intercept, shape_slope, USE_VIENNA, maxdist, cmd_pk, DMS_file )
% run_rna_structure_with_bootstrap( NUM_BOOTSTRAP, seq_file, bps, offset, nres, EX_file, SHAPE_file, temperature, experimental_offset, zscore_scaling, shape_intercept, shape_slope, USE_VIENNA, maxdist );
%
% make sure to set your path in get_exe_dir.m
%

if ~exist( 'maxdist','var' ); maxdist = 0; end;
if ~exist('cmd_pk','var'); cmd_pk = 0; end;
if ~exist('DMS_file','var'); DMS_file = ''; end;


[structure, bpp, ct_file, command ] = run_rna_structure_with_EX_and_SHAPE( seq_file, temperature, experimental_offset, zscore_scaling, EX_file, SHAPE_file, shape_intercept, shape_slope, USE_VIENNA, maxdist, cmd_pk, DMS_file );
system( ['rm ',ct_file ] );

if NUM_BOOTSTRAP == 0; return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boostrap loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_bpp = zeros( nres, nres, NUM_BOOTSTRAP );
for i = 1:NUM_BOOTSTRAP;   structure_boot{i} = ''; end

% this is set up to work with the Matlab parallelization toolbox.
if exist( 'matlabpool' ) && parallelization_exists()
    if matlabpool( 'size' ) == 0;    res = findResource;  matlabpool( res.ClusterSize ) ; end
    parfor n = 1:NUM_BOOTSTRAP
        [structure_boot{n}, all_bpp(:,:,n) ] = main_loop( n, seq_file, temperature, experimental_offset, zscore_scaling, EX_file, SHAPE_file, shape_intercept, shape_slope, USE_VIENNA, maxdist, cmd_pk, DMS_file );
    end
else
    for n = 1:NUM_BOOTSTRAP
        [structure_boot{n}, all_bpp(:,:,n) ] = main_loop( n, seq_file, temperature, experimental_offset, zscore_scaling, EX_file, SHAPE_file, shape_intercept, shape_slope, USE_VIENNA, maxdist, cmd_pk, DMS_file );
    end
end

char( structure )
char( structure_boot )
bpp = squeeze( sum( all_bpp, 3 ) ) / NUM_BOOTSTRAP;


% Plot will be made outside
%make_plot( bpp, nres, bps, offset );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [structure_boot, bpp_boot ] = main_loop( n, seq_file, temperature, experimental_offset, zscore_scaling, EX_file, SHAPE_file, shape_intercept, shape_slope, USE_VIENNA, maxdist, cmd_pk, DMS_file )

EX_file_boot = '';
SHAPE_file_boot = '';
DMS_file_boot = '';

if ~isempty(EX_file);
    EX_file_boot = ['EX_',num2str(n),'.txt'] ;
    create_bootstrap_EX( EX_file, EX_file_boot );
end
if ~isempty(SHAPE_file);
    SHAPE_file_boot = ['SHAPE_',num2str(n),'.txt'];
    create_bootstrap_SHAPE( SHAPE_file, SHAPE_file_boot );
end
if ~isempty(DMS_file);
    DMS_file_boot = ['DMS_',num2str(n),'.txt'];
    create_bootstrap_SHAPE( DMS_file, DMS_file_boot );
end

[ structure_boot, bpp_boot, ct_file ] = run_rna_structure_with_EX_and_SHAPE( seq_file, temperature, experimental_offset, zscore_scaling, EX_file_boot, SHAPE_file_boot, shape_intercept, shape_slope, USE_VIENNA, maxdist, cmd_pk, DMS_file_boot );

if ~isempty(EX_file_boot)    & exist( EX_file_boot, 'file' );   delete( EX_file_boot ); end
if ~isempty(SHAPE_file_boot) & exist( EX_file_boot, 'file' );   delete( SHAPE_file_boot ); end
if ~isempty(DMS_file_boot)   & exist( DMS_file_boot, 'file');   delete( DMS_file_boot ); end
delete( ct_file );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_plot( bpp, nres, bps, offset )

image( 128 * bpp );
colormap( 1 - copper( 128 ) );
hold on

for i = 1:size(bps,1)
    rectangle( 'Position', [bps(i,1)-0.5, bps(i,2)-0.5, 1.0, 1.0], 'edgecolor','r');
end
plot( [1:nres],[1:nres], 'k' );

seqpos = offset + [1 : nres ];
hold off

gps = find( mod( seqpos, 10 ) == 0 );
set(gca,'xtick',gps,'xticklabel',seqpos(gps ),'xgrid','on');
set(gca,'ytick',gps,'yticklabel',seqpos(gps ),'ygrid','on');
set(gca,'ydir','normal');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  create_bootstrap_EX( EX_file, EX_file_boot )

ex = load( EX_file );
nres = length(ex);
boot_idx = randi( nres, 1, nres );

for i = 1:nres
    ncontrib = sum( boot_idx == i  );
    ex_boot(:,i) = ex(:,i) * ncontrib;
end

save( EX_file_boot, 'ex_boot', '-ascii');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  create_bootstrap_SHAPE( SHAPE_file, SHAPE_file_boot )

shape = load( SHAPE_file );
nres = length(shape);
boot_idx = randi( nres, 1, nres );

boot_idx = sort( boot_idx );
shape_boot = shape( boot_idx, :);

fid = fopen( SHAPE_file_boot, 'w' );
for k = 1:nres
    fprintf( fid, '%d %9.4f\n', shape_boot(k,1), shape_boot(k,2) );
end
fclose( fid );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [structure, bpp, ct_file, command ] = run_rna_structure_with_EX_and_SHAPE( seq_file, temperature, experimental_offset, zscore_scaling, EX_file, SHAPE_file, shape_intercept, shape_slope, USE_VIENNA, maxdist, cmd_pk, DMS_file );

if USE_VIENNA
    [structure, bpp, ct_file, command  ] = run_vienna_with_EX_and_SHAPE( seq_file, temperature, experimental_offset, zscore_scaling, EX_file, SHAPE_file, shape_intercept, shape_slope, USE_VIENNA );
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic command
if cmd_pk == 0;
    cmd = 'Fold';
else
    cmd = 'ShapeKnots';
end;
EXE = [get_exe_dir(), cmd];


if ~isempty( EX_file )
    ct_file = [ EX_file(end-6:end),'.temp.ct'];
elseif ~isempty( SHAPE_file )
    ct_file = [ SHAPE_file(end-6:end),'.temp.ct'];
else
    ct_file = [ seq_file(end-6:end),'.temp.ct'];
end

if cmd_pk == 0;
    command = [EXE,' ',seq_file,' ',ct_file,' -T ',num2str( 273.15 + temperature, '%6.2f' ) ];
else
    command = [EXE,' ',seq_file,' ',ct_file ];
    
end;

if ( maxdist > 0 ); command = [ command, ' -md ', num2str( maxdist) ];  end;

if ~isempty( EX_file ); command = [ command, ' -x ', EX_file ]; end;
if ~isempty( SHAPE_file ); command = [ command, ' -sh ', SHAPE_file ]; end;
if ~isempty( DMS_file ); command = [ command, ' -dms ', DMS_file ]; end;

if  ~isempty( shape_intercept); command = [ command, ' -si ', num2str( shape_intercept ) ]; end;
if  ~isempty( shape_slope ); command = [ command, ' -sm ', num2str( shape_slope ) ]; end;

if experimental_offset ~= 0.0
    command = [ command, ' -xo ', num2str( experimental_offset ) ];
end

if zscore_scaling ~= 1.0
    command = [ command, ' -xs ', num2str( zscore_scaling ) ];
end


fprintf( 'COMMAND:   %s\n',command )
system( command );
[structure, bpp] = load_ct_file( ct_file );

%system( ['rm ',ct_file ] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [structure, bpp, pfs_file, command  ] = run_vienna_with_EX_and_SHAPE( seq_file, temperature, experimental_offset, zscore_scaling, EX_file, SHAPE_file, shape_intercept, shape_slope, USE_VIENNA )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic command
EXE = 'python ../../../external/ViennaRNA-1.8.4/Utils/viennafold.py';

if ~isempty( EX_file );
    pfs_file = [ EX_file(end-6:end),'.temp.pfs'];
elseif ~isempty( SHAPE_file )
    pfs_file = [ SHAPE_file(end-6:end),'.temp.pfs'];
else
    pfs_file = [ seq_file(end-6:end),'.temp.pfs'];
end;

command = [EXE,' ',seq_file,' ',pfs_file];% currently no temperature!!! ,' -T ',num2str( 273.15 + temperature, '%6.2f' ) ];

if temperature ~= 24.0; fprintf( 'temperature disabled for vienna!\n' );end;

if ~isempty( EX_file ); command = [ command, ' -x ', EX_file ]; end;
if ~isempty( SHAPE_file ); command = [ command, ' --sh ', SHAPE_file ]; end;

if  ~isempty( shape_intercept); command = [ command, ' --si ', num2str( shape_intercept ) ]; end;
if  ~isempty( shape_slope ); command = [ command, ' --sm ', num2str( shape_slope ) ]; end;

if experimental_offset ~= 0.0;  fprintf( 'experimental offset disabled for vienna!\n' ); end;

if zscore_scaling ~= 1.0;  fprintf( 'experimental offset disabled for vienna!\n' ); end;


fprintf( 'COMMAND:   %s\n', command )
[error, output ] = system( command );

fprintf( output )

[t,r] = strtok( output, '\n' );
[t,r] = strtok( r, ':' );
[t,r] = strtok( r, ',' );

structure = t(3:end);
bpp = load( pfs_file );

%clf;
%image( 100*bpp )
%pause;

%[structure, bpp] = load_pfs_file( pfs_file );

