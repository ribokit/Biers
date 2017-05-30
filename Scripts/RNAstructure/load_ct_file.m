function [structure,bpp] = load_ct_file( ct_file )
% [structure,bpp] = load_ct_file( ct_file )
% only loads first structure!!
%
%
% (C) Das lab, Stanford University 2011-2015

fid = fopen( ct_file );

line = fgetl( fid );
nres = str2num( strtok( line, ' ' ) );
bpp = zeros( nres,nres );

structure = '';
partner_prv = 0;
helix_l = 0;
helix_map = [];
helix_map_ct = 0;

for count = 1:nres;
    line = fgetl( fid );
    
    for j = 1:5;  [t,line] = strtok( line, ' ' ); end;
    partner = str2num( t );
    
    if abs(partner_prv - partner) == 1 && partner ~= 0;
        helix_l = helix_l + 1;
    elseif helix_l ~= 0 && partner_prv ~= 0;
        if partner_prv >= count-1;
            helix_map_ct = helix_map_ct + 1;
            helix_map(helix_map_ct,:) = [count-1-helix_l,partner_prv+helix_l,helix_l+1];
        end;
        helix_l = 0;
    end;
    
    if ( partner == 0 );
        symbol = '.';
    else
        bpp( count, partner ) = 1.0;
        if ( partner > count );
            symbol = '(';
        else
            symbol = ')';
        end;
    end;
    structure = [ structure, symbol ];
    
    partner_prv = partner;
    
end;


fclose( fid );

clash_score = zeros(1,size(helix_map,1));
for i = 1:size(helix_map,1)
    for j = 1:size(helix_map,1)
        if (helix_map(i,1)<helix_map(j,1) && helix_map(i,2)<helix_map(j,2) && helix_map(i,2) > helix_map(j,1)) ...
                || (helix_map(i,1)>helix_map(j,1) && helix_map(i,2)>helix_map(j,2) && helix_map(j,2) > helix_map(i,1))
            clash_score(i) = clash_score(i)+1;
        end;
    end;
end;


clash_score_flag = 0;
for i = 1:length(clash_score)
    if clash_score(i) > 1;
        structure = pk_braket_substitute(structure, helix_map(i,:));
        clash_score_flag = 1;
    end;
end;
if ~clash_score_flag && any(clash_score ~= 0);
    clash_odd =[];
    clash_odd_flag = 0;
    for i = 1:length(clash_score)
        if clash_score(i) == 1;
            if clash_odd_flag == 0;
                clash_odd_flag = 1;
                clash_odd = [clash_odd, i];
            else
                clash_odd_flag = 0;
            end;
        end;
    end;
    for i = 1:length(clash_odd)
        structure = pk_braket_substitute(structure, helix_map(clash_odd(i),:));
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str_return = pk_braket_substitute(str_input, helix_map_sub)
str_return = [str_input(1:(helix_map_sub(1)-1)), repmat('[',1,helix_map_sub(3)), ...
    str_input((helix_map_sub(1)+helix_map_sub(3)):(helix_map_sub(2)-helix_map_sub(3))),...
    repmat(']',1,helix_map_sub(3)), str_input((helix_map_sub(2)+1):end)];

