function EXE_DIR = get_RNAstructure_exe_dir()

% make sure to change EXE_DIR to your local path
% the following line is just a dummy path
% it is used to locate the RNAstructure executables, e.g. Fold and ShapeKnot

% use full path, e.g. "/User/yourname/Desktop/RNAstructure/exe/", with trailing slash

rna_structure_datapath = getenv( 'DATAPATH' );
if length( rna_structure_datapath ) == 0 | ~exist( rna_structure_datapath, 'dir' )
    error( sprintf( '\n\nYou must install RNAstructure and set the DATAPATH variable in your .bashrc.\nIf you did that, make sure you start matlab from Terminal with %s/bin/matlab\n', matlabroot ) );
end

[EXE_DIR,datapathname] = fileparts(rna_structure_datapath);

% user may have included a trailing slash in 'data_tables/':
if length( datapathname ) == 0; [ EXE_DIR, datapathname ] = fileparts( EXE_DIR ); end;

assert( strcmp( datapathname, 'data_tables' )  )
EXE_DIR = [EXE_DIR, '/exe/'] ;
if ~exist( EXE_DIR, 'dir' ) open( mfilename ); return; end;
