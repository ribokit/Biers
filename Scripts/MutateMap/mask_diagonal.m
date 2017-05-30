function a_mask = mask_diagonal( a, mask_diag);
% a_mask = mask_diagonal( a, mask_diag);
%
% INPUTS
%  a        = square matrix
% mask_diag = remove points with i-j within this number. (-1 means no mask)
%
% (C) R. Das, Stanford University, 2017

if ( nargin< 2 ) help( mfilename ); return; end

a_mask = a;
if ( mask_diag >= 0 )
    [idx_i, idx_j ] = ndgrid( [1:size(a,1)], [1:size(a,2)] );
    diag_pts = find( abs( 1 + idx_i - idx_j ) <= mask_diag ); 
    a_mask( diag_pts ) = 0;
 end
