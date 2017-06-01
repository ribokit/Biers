function [Zfinal,structure] = m2net( Z, sequence, structure, offset )
% [Zfinal,structure] = m2net( r )
% [Zfinal,structure] = m2net( Z, r )
% [Zfinal,structure] = m2net( Z, sequence, structure, offset )
%
% Extremely simple convolutional neural net to identify
%  helices in mutate-and-map (M2) data.
%
% INPUTS
%
%  r = RDAT Object
%  Z = Z-score matrix
%  sequence  = sequence (use lowercase for flanking regions to ignore)
%  structure = secondary structure in dot-parens rotation for comparison
%  offset    = offset to add to 1,2,... to get conventional numbering.
%
%
% Output
% Zfinal = final layer of net -- scored base-pairs. 
%
% (C) R. Das, Stanford University, 2017

% options setting [could spin this out into a user-inputtable option
opts.mask_diagonal = 7;
opts.mask_flanking_lowercase_sequence = true;
opts.cross_diagonal_filter_stride = 5;
opts.cross_diagonal_relu_bias = 1;
opts.minimum_helix_length = 3;
opts.filter_for_single_partner = true;

if ( nargin < 1 ) help( mfilename ); end;
if nargin == 1 & isobject( Z )
    r = Z;  % first arg
    Z = output_Zscore_from_rdat( [], r, [], [], 1, 1 );
    sequence = r; 
end
if isobject( sequence )
    r = sequence;
    sequence = r.sequence;
    structure = r.structure;
    offset = r.offset;
end
L = length( sequence );

% prep one-hot arrays (tiled vertical and horizontal) of isA, isC, isG, isU
rnachars = 'ACGU';
for i = 1:length( rnachars )
    ischar = sequence == rnachars(i);
    [Dx{i},Dy{i}] = ndgrid( ischar, ischar );
end

% figure out where Watson-Crick pairs might occur.
base_pairs = {'AU','CG','GU'};
for i = 1:length( base_pairs ); base_pairs = [base_pairs, fliplr(base_pairs{i})]; end

bp = zeros( L, L );
for i = 1:length( base_pairs )
    base_pair = base_pairs{i};
    base_pair = [strfind(rnachars,base_pair(1)), strfind(rnachars,base_pair(2))];
    bp = bp + Dx{base_pair(1)} .* Dy{base_pair(2)};
end
%image( bp*30 );

% figure out helices of length n = 1 ... 10.
N = 10;
% easiest to first mark ends of helices
ends_helix{1} = bp;
for n = 2:N
    ends_helix{n} = ends_helix{n-1} .* circshift(bp,[1 -1]);
end
starts_helix{1} = bp;
for n = 2:N
    starts_helix{n} = starts_helix{n-1} .* circshift(bp,[1 -1]);
end
% then mark each helix (propagating back from the end)
for n = 1:N
    in_helix{n} = ends_helix{n};
    for k = 1:(n-1)
        in_helix{n} = max( in_helix{n}, circshift(ends_helix{n},[-k,k]) );
    end
end

% mask more of the diagonal
Zmask = mask_diagonal( Z, opts.mask_diagonal );
if ( opts.mask_flanking_lowercase_sequence )
    % also mask out 'flanking' sequences, given in lowercase
    flankpos = find( lower(sequence) == sequence );
    Zmask( flankpos, : ) = 0;
    Zmask( :, flankpos ) = 0;
end
show_2dmap( 50*-Zmask, structure, offset );

% apply 5x5 filter to look for stripe.
S = opts.cross_diagonal_filter_stride; % filter size
B_bias = 0; %-0.1;
B = ( diag(ones(S,1) ) + B_bias)/S; % filter
B = fliplr( B ); % cross-diagonal
Z2 = filter2( B, -Zmask );
Z3 = max( (Z2 + Z2') - opts.cross_diagonal_relu_bias, 0 );
show_2dmap( 50*Z3, structure, offset ); 

% now filter for helices > 2bp.
Z4 = Z3.*tril(ends_helix{ opts.minimum_helix_length }); 

% remove singlets
Z5 = Z4;
Z5( ~circshift(Z4,[-1, 1]) & ~circshift(Z4,[1, -1]) ) = 0;

% filter so that no position has more than one partner.
[m,idx] = sort( Z5(:) );
Z6 = 0*Z5;
used = zeros( 1, L );
for q = fliplr(idx')
    [i,j] = ind2sub( size(Z5), q );
    if ( ~used( i ) & ~used( j ) )
        Z6(i,j) = Z5(i,j); used(i) = 1; used(j) = 1;
    end
end

% filter singlets again
Z7 = Z6;
Z7( ~circshift(Z7,[-1, 1]) & ~circshift(Z7,[1, -1]) ) = 0;


%subplot(1,2,1);
show_2dmap( 200*Z7, structure, offset );

%hold on; plot( sum(how_competitive) ); plot( sum(how_competitive') ); hold off;
% apply 'smooth' filter to identify less punctate helices.
%Y = -Zmask;
%Y2 = smooth2d(Y,5);
%subplot(1,2,2);
%show_2dmap( 50*Y2, structure, offset );

Zfinal = Z7;

% convert to structure in dot-parens notation
[ix,jx] = ind2sub( size( Zfinal ), find( Zfinal' > 0 ) );
bps = [ix,jx];
structure = convert_bps_to_structure( bps, length(Zfinal) );
