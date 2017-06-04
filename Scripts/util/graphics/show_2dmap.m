function show_2dmap( Z, structure, offset, maxZ, show_colorbar );
% show_2dmap( Z, structure, offset );
%
%   Makes plot of square matrix of scores, e.g., 
%       -Z scores, base pair probabilities,
%       or bootstrapping probabilities.
%   Pulled out of cluster_z_scores
%
% INPUTS:
% Z         = NxN square matrix of *positive* counts (you may need to 
%               negate a Z-score matrix)
% structure = structure string in dot parens notation -- will show up 
%               as squares.
% offset    = value to add to 1, 2, ... N to get conventional numbering
% maxZ      = [optional] maximum value of Z to show (default is [], no
%              scaling)
% show_colorbar = show colorbar to side (default 0).
%
if nargin < 1; help( mfilename ); return; end;
if ( size( Z, 1 ) ~= size( Z, 2 ) ); 
    if size( Z, 1 ) == size( Z, 2 ) + 1;
        Z = Z(:,2:end);
        fprintf( 'Warning: assuming you have an M2 file, but need a square plot. Not showing wild type!' );
    else
        fprintf( 'Z needs to be square' );
        return;
    end
end
if ~exist( 'maxZ' ) || isempty( maxZ );
    image( Z' );
    maxZ = 100;
else
    if ( maxZ < 0 ); Z = -Z; maxZ = -maxZ; end;
    imagesc( Z', [0, maxZ] );
end

colormap(gca,  1 - gray(100) );
set_axes( Z, offset );
show_bps( structure );
if show_colorbar
    h = colorbar;
    set( h,'ticks',[0, maxZ]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_bps( structure, USE_DOTS );

if length( structure ) == 0; return; end;
if ~exist( 'USE_DOTS' ); USE_DOTS = 0; end;

hold on;
bps = convert_structure_to_bps( structure );
bps = [ bps; bps(:,[2 1] ) ];
if USE_DOTS
  plot( bps(:,1), bps(:,2), 'ko','markersize',2,'markerfacecolor','w','linewidth',1 );
else
  for i = 1:size( bps, 1 )
    h = rectangle( 'Position', [ bps(i,1)-0.5, bps(i,2)-0.5, 1, 1], ...
		   'edgecolor', [1 0.5 0.5],...
		   'linewidth',0.5 );  
  end
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_axes( d, offset, axiscolor )
if ~exist( 'axiscolor' );  axiscolor = 'k'; end;
nres = length( d );
seqpos = [1:nres] + offset;
gp = find( mod(seqpos,20) == 0 );
set( gca,'xtick', gp, 'xticklabel',  seqpos( gp ), 'ytick', gp, 'yticklabel', seqpos( gp ) );
set( gca,'xgrid','on','ygrid','on','fontsize',12,'fontweight','bold','xcolor',axiscolor,'ycolor',axiscolor );
xticklabel_rotate;
hold on
plot( [1:nres],[1:nres],axiscolor);
hold off
set(gcf,'color','white');
xlabel( 'Mapped Position' );
ylabel( 'Mutation Position' );
axis image;
