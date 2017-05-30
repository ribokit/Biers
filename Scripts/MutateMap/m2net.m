function m2net( Z, sequence, structure, offset )
%
% m2net( Z, r )
% m2net( Z, sequence, structure, offset )
%
% Extremely simple convolutional neural net to identify
%  helices in mutate-and-map (M2) data.
%
%


if nargin == 2 & isobject( sequence )
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
Zmask = mask_diagonal( Z, 7 );

% apply 5x5 filter to look for stripe.
S = 7; % filter size
B = diag( ones(S,1) )/S; % filter
B = fliplr( B ); % cross-diagonal
Z2 = filter2( B, -Zmask );
Z3 = max( (Z2 + Z2')-1, 0 );

% now filter for helices > 2bp.
Z4 = Z3.*tril(ends_helix{3}); 

% remove singlets
Z5 = Z4;
Z5( ~circshift(Z4,[-1, 1]) & ~circshift(Z4,[1, -1]) ) = 0;

%subplot(1,2,1);
show_2dmap( 20*Z5, structure, offset );

% apply 'smooth' filter to identify less punctate helices.
Y = -Zmask;
Y2 = smooth2d(Y,5);
%subplot(1,2,2);
%show_2dmap( 50*Y2, structure, offset );

