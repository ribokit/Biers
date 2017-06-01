function [structure,bpp] = load_ct_file( ct_file )
% [structure,bpp] = load_ct_file( ct_file )
% only loads first structure!
%
% INPUTS 
%  ct_file = file in ct format, as is used by, e.g., the RNAstructure
%  
%
% (C) Das lab, Stanford University 2011-2015, 2017


% initialize structure with '.'
structure = '';
for i = 1:nres; structure = [structure,'.']; end;
left_delim  = '([{abcdefghijklmnopqrstuvwxyz';
right_delim = ')]}abcdefghijklmnopqrstuvwxyz';
if ( size( helix_map, 1 ) == 0 ) return; end;

% read ct_file
fid = fopen( ct_file );

line = fgetl( fid );
nres = str2num( strtok( line, ' ' ) );
bpp = zeros( nres,nres );

helix_map = [];

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

stems = parse_stems_from_bps( bps )
for i = 1:length(stems)
    stems{i}
    helix_map = [helix_map; stems{i}(1,1), stems{i}(1,2), size( stems{i}, 1 ) ];
end


    
% order helices by length
[~,idx] = sort( helix_map(:,3) );
idx = fliplr( idx );

% "layers" are helices that can be connected by (), then by [], then by {},
% then by aa, ...
helix_layers = {};
for n = 1:length( idx )
    i = idx( n );
    found_layer = 0;
    test_helix = helix_map(i,:);
    for j = 1:min(3,length(helix_layers))
        if ( not_pseudoknotted( helix_layers{j}, test_helix ) )
            helix_layers{j} = [ helix_layers{j}, test_helix ];
            found_layer = 1;
        end
    end
    if ( ~found_layer ) % new layer
        helix_layers = [helix_layers, helix_map(i,:) ];
    end
end

for j = 1:length( helix_layers )
    helix_layer = helix_layers{j};
    for i = 1:size( helix_layer, 1 )
        structure = pk_bracket_substitute(structure, helix_layer(i,:), left_delim(j), right_delim(j) );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = not_pseudoknotted( helix_layer, test_helix )

ok = 1;
for i = 1:size( helix_layer, 1 )
    % "crossing" base pairs
    if ( helix_layer(i,1) < test_helix(1) & ...
         helix_layer(i,2) < test_helix(2) & ...
         helix_layer(i,2) > test_helix(1) ) ...
         | (helix_layer(i,1) > test_helix(1) & ...
            helix_layer(i,2) > test_helix(2) & ...
            test_helix(2) > helix_layer(i,1) )
        ok = 0;
        return;
    end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str_return = pk_bracket_substitute(str_input, helix_layer_sub, left_delim, right_delim)
str_return = [...
    str_input(1:(helix_layer_sub(1)-1)), ...
    repmat(left_delim,1,helix_layer_sub(3)), ...
    str_input((helix_layer_sub(1)+helix_layer_sub(3)):(helix_layer_sub(2)-helix_layer_sub(3))),...
    repmat(right_delim,1,helix_layer_sub(3)), ...
    str_input((helix_layer_sub(2)+1):end)];

