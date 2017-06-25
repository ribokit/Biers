function [VARNA_DIR, VARNA_JAR] = get_varna()
%
% Supplies local path, e.g. "/Users/yourname/src/VARNA.jar", based on
% enironment variable $VARNA.
%
% Originally this gave a url like
% "https://rmdb.stanford.edu/site_data/VARNA.jar", but
% no web browser supports Java anymore.
%
% (C) R. Das, Stanford University, 2017

varna = getenv( 'VARNA' );
if length( varna ) == 0 | ( ~exist( varna, 'file' ) & isempty(strfind( varna, 'http' )) )
    error( sprintf( '\n\nYou must download VARNA.jar and include a line like "export VARNA=/path/to/VARNA.jar" in your .bashrc.\nIf you did that, make sure you start matlab from Terminal with %s/bin/matlab\n', matlabroot ) );
end

[VARNA_DIR, VARNA_JAR, EXT] = fileparts( varna );
VARNA_DIR = [VARNA_DIR,'/'];
VARNA_JAR = [VARNA_JAR, EXT];