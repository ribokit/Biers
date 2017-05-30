function show_2dmap( Z, structure, offset );
%
% Pulled out of cluster_z_scores
%
%

image( Z' );
colormap( 1 - gray(100) );
set_axes( Z, offset );
show_bps( structure );

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
set( gca,'xticklabelrotation', 90 );
hold on
plot( [1:nres],[1:nres],axiscolor);
hold off
set(gcf,'color','white');
xlabel( 'Sequence Position' );
ylabel( 'Mutation Position' );
axis image;
