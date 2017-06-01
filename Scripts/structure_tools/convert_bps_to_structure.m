function structure = convert_bps_to_structure( bps, nres )
% structure = convert_bps_to_structure( bps, nres )
%
%  Note: wraps around load_ct_file()
%
% INPUTS
%  bps  = list of base pairs (i,j) (Nbp X 2 matrix )
%
%  (C) R. Das, Stanford University, 2009-2015, 2017

if ( nargin < 2 ) help( mfilename ); return; end;

if length( bps ) > 0 & size( bps, 2 ) ~= 2
    fprintf( 'bps must be Nbps x 2 array\n' );
    return;
end

bps_all = zeros(nres,2);
for i = 1:nres
    bps_all(i,2) = i;
    if find(bps(:,2) == i) ~= 0;
        bps_all(i,1) = bps(find(bps(:,2) == i),1);
    end;
    if find(bps(:,1) == i) ~= 0;
        bps_all(i,1) = bps(find(bps(:,1) == i),2);
    end;
end;

fid = fopen('bps_temp.txt','w');
fprintf(fid,'  %d  ENERGY = 0  temp\n',nres);
for i = 1:nres
    fprintf(fid,'    %d A       %d    %d    %d    %d\n',i,i-1,i+1,bps_all(i,1),bps_all(i,2));
end;
fclose(fid);

structure = load_ct_file( 'bps_temp.txt' );
system('rm bps_temp.txt');
    