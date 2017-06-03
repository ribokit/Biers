function domains = define_domains( stems, loop2_cutoff )
% domains = define_domains( stems, loop2_cutoff );
%
% Groups stems into domains -- stems that
%   are connected by 2-way junctions (internal loops) with
%   total length less than L are grouped together.
%
% INPUTS:
% stems = cell of Nx2 arrays of base pairs (use, e.g., parse_stems or
%           parse_stems_from_bps )
% loop2_cutoff = after this value, stems are considered
%          in different domains. loop2_cutoff is applied to sum of the lengths
%          two single-stranded segments connecting helices [default 8]
%
% OUTPUTS:
% domains = assignment of each stem to a domain
%
if nargin < 1; help( mfilename ); return; end;
if ~exist( 'loop2_cutoff', 'var' ) loop2_cutoff = 8; end;

% set up connection graph.
N = length(stems);
connections = diag( ones(1,N ) );
for i = 1:N
    for j = 1:i-1
        connections( i, j ) = check_connected( stems, i, j, loop2_cutoff );
        connections( j, i ) = connections( i, j );
    end
end

% find connected components of this graph:
domains = zeros( 1, N );
for n = 1:N
    if ( domains( n ) ~= 0 ) continue; end;
    partners = find( connections(n,:) > 0);
    domain = max( domains( partners ) );
    if ( domain == 0 ) domain = max( domains ) + 1; end;
    domains( partners ) = domain;
end

function connected = check_connected( stems, i, j, loop2_cutoff );
% should be nested:  (((((..(((((.....)))))...)))))
%                    IIIII  JJJJJ     JJJJJ   IIIII
%                    11111  11111     22222   22222
%                    1 end  1 end     end 1   end 1
%                        L_ij              L_ji
%
% Note: this assumes that stems are different stems drawn from a
%  legitimate secondary structure -- i.e., the stems
%  cannot share any base pairs.
%
if ( i == j ); connected = 1; return; end;
stemI = stems{i};
stemJ = stems{j};
% order the stems so that 5' strand of B is after 5' strand of A
if  stemJ(1,1) < stemI(1,1)
    connected = check_connected( stems, j, i, loop2_cutoff );
    return;
end
assert( stemJ(1,1) >= stemI(end,1) );
% Nested
if ( stemI(1,2) < stemJ(1,2) );
    connected = 0; return
end;
assert( stemI(end,2) >= stemJ(1,2) );
% check for any intervening stems
for n = 1:length(stems)
    if ( n == i | n == j ); continue; end;
    if ( stems{n}(1,1) > stemI(end,1) & stems{n}(end,1) < stemJ(end,1) )
        connected = 0; return;
    end
    if ( stems{n}(1,1) > stemJ(1,2) & stems{n}(end,1) < stemI(1,2) )
        connected = 0; return;
    end
end

L_ij = ( stemJ(1,1) - stemI(end,1) - 1);
L_ji = ( stemI(end,2) - stemJ(1,2) - 1 );
L = L_ij + L_ji;
connected = ( L <= loop2_cutoff );
return;




