function DMS_profile = get_DMS_profile( r, r_nomod );
% DMS_profile = get_DMS_profile( r, r_nomod );
% DMS_profile = get_DMS_profile( {r, r_nomod} );
%
% DMS profile with NaN's at G, U.
%
% NOTE: Currently only really works with M2seq-DMS rdats,
%  whose first line is the 1D DMS profile, and
%  which has values for all sequence positions.
%
% But could be easily generalized.
%
% (C) R. Das, Stanford University, 2017

if iscell( r )
    r_nomod = r{2};
    r = r{1};
end
if ~isobject( r ) & ischar( r )
    r = read_rdat_file( r );
end
if exist( 'r_nomod', 'var' ) & ~isobject( r_nomod ) & ischar( r_nomod ) & length( r_nomod ) > 0
    r_nomod = read_rdat_file( r_nomod );
end

DMS_profile = r.reactivity(:,1);
DMS_pos = [strfind(r.sequence,'A'),strfind(r.sequence,'C')]; %,strfind(r.sequence,'a'),strfind(r.sequence,'c')];
non_DMS_pos = setdiff( [1:size( r.reactivity, 1 )], DMS_pos );

if exist( 'r_nomod', 'var' ) & isobject( r_nomod )
    DMS_profile = DMS_profile - r_nomod.reactivity(:,1);
else
    % do background subtraction based on closest non-DMS pos
    for i = 1:length( DMS_pos )
        [~,closest_idx] = min( abs(i - non_DMS_pos) );
        DMS_profile( i ) = DMS_profile( i ) - DMS_profile( non_DMS_pos( closest_idx ) );
    end
end

DMS_profile( non_DMS_pos ) = NaN;
