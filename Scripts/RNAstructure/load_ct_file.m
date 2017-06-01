function [structure,bpp] = load_ct_file( ct_file )
% [structure,bpp] = load_ct_file( ct_file )
% only loads first structure!
%
% INPUTS 
%  ct_file = file in ct format, as is used by, e.g., the RNAstructure
%  
%
% (C) Das lab, Stanford University 2011-2015, 2017


% read ct_file
fid = fopen( ct_file );

line = fgetl( fid );
nres = str2num( strtok( line, ' ' ) );
bpp = zeros( nres,nres );

% initial fill-in of helix_map whose rows are:
%  begin-helix, end-helix, length of helix
%  where begin-helix and end-helix are base paired.
bps = [];
for count = 1:nres;
    line = fgetl( fid );
        for j = 1:5;  [t,line] = strtok( line, ' ' ); end;
    partner = str2num( t );
    if partner ~= 0 & partner > count
        bps = [bps;count, partner ];
    end
end;
fclose( fid );

structure = convert_bps_to_structre( bps, nres );
