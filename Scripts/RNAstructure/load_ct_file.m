function [structure,bpp,bps] = load_ct_file( ct_file )
% [structure,bpp, bps] = load_ct_file( ct_file )
% only loads first structure!
%
% INPUTS 
%  ct_file = file in ct format, as is used by, e.g., the RNAstructure
%
% OUTPUTS
%  structure = (string of length Nres) dot bracket notation, same length as RNA sequence
%  bpp       = [Nres x Nres] base pair matrix (0 and 1)
%  bps       = [Nbp x 2] list of base pairs
%  
% (C) Das lab, Stanford University 2011-2015, 2017, 2024


% read ct_file
fid = fopen( ct_file );

line = fgetl( fid );
nres = str2num( strtok( line ) );
bpp = zeros( nres,nres );

% initial fill-in of helix_map whose rows are:
%  begin-helix, end-helix, length of helix
%  where begin-helix and end-helix are base paired.
bps = [];
for count = 1:nres;
    line = fgetl( fid );
    for j = 1:5;  [t,line] = strtok( line ); end;
    res = [];
    if length(line) > 0; res = str2num( strtok( line ) ); end;
    partner = str2num( t );
    if isempty(res); res = count; end;
    if partner ~= 0 & partner > res
        bps = [bps;res, partner ];
    end
end;
fclose( fid );

structure = convert_bps_to_structure( bps, nres );

for i = 1:size( bps, 1); 
    bpp( bps(i,1), bps(i,2) ) = 1; 
    bpp( bps(i,2), bps(i,1) ) = 1; 
end
    
