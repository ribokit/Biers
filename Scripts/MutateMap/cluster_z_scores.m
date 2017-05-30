function  cluster_z_scores( zscore_matrix, structure, offset, print_mode  ); 
% cluster_z_scores( zscore_matrix, structure, offset, print_mode )
%
% INPUTS
%  zscore_matrix = [REQUIRED] square Z matrix (or filename with that matrix )
%  structure     = structure in dot-bracket notation
%  offset        = offset to add to position to get conventional numbering
%  print_mode    = 1: show plots in fig(1),fig(2) [default] 
%                  0: show final clusters but do not go to a new figure panel.
%
% (C) R. Das, 2011-2013, 2017

if nargin == 0;  help( mfilename ); return; end;

if ischar( zscore_matrix ); 
  d = load( zscore_matrix );
else  
  d = zscore_matrix;
end

if ~exist( 'offset', 'var' ) offset = 0; end;
if ~exist( 'structure', 'var' ); structure = ''; end;
if ~exist( 'print_mode', 'var' ) print_mode = 1; end;
d = abs( smooth2d(d,1) );

if ( print_mode )
    figure(1)
    show_2dmap( d*20, [], offset );
    figure(2)
end

NRES = length( d ) ;
cluster_assignment = zeros( NRES, NRES );

count = 0;
CUTOFF = 1;
for i = 1:NRES
  for j = 1:NRES
    if d(i,j) > CUTOFF
      count = count + 1;
      cluster_assignment(i,j) = count;
    end
  end
end


NPASS = 5;
for n = 1:NPASS
    for i = 1:NRES
        for j = 1:NRES
            if ( cluster_assignment( i, j ) )
                cluster_assignment = check_neighbor( cluster_assignment, i, j, i-1, j );
                cluster_assignment = check_neighbor( cluster_assignment, i, j, i+1, j );
                cluster_assignment = check_neighbor( cluster_assignment, i, j, i,   j-1);
                cluster_assignment = check_neighbor( cluster_assignment, i, j, i,   j+1 );
                
                %cluster_assignment = check_neighbor( cluster_assignment, i, j, i-1, j+1 );
                %cluster_assignment = check_neighbor( cluster_assignment, i, j, i+1, j-1 );
                
                cluster_assignment = check_neighbor( cluster_assignment, i, j, j,   i );
            end
        end
    end
    

end

if ( print_mode )
    figure(1)
    subplot(1,2,1)
    show_2dmap( cluster_assignment, structure, offset );
    cmap = jet( count );
    cmap = cmap( randperm( count ), : );
    cmap( 1,: ) = 0.0;
    colormap( cmap );
end

% filter out weak clusters.
CLUSTER_CUTOFF = 8;
NUM_MUT_CUTOFF = 3; % disabled
ZSCORE_SUM_CUTOFF = 0.0; % disabled
MAX_ZSCORE_CUTOFF = 0.0; % disabled
MEAN_ZSCORE_CUTOFF = 0.0; % disabled
AVG_HITS_PER_MUT_CUTOFF = 1000; % disabled
CLUSTERS_PER_MUT_CUTOFF = 5; % important
FORCE_SYMM = 1;

cluster_assignment_plot = get_cluster_filter( d, cluster_assignment, CLUSTER_CUTOFF, ZSCORE_SUM_CUTOFF,MAX_ZSCORE_CUTOFF, MEAN_ZSCORE_CUTOFF, NUM_MUT_CUTOFF, AVG_HITS_PER_MUT_CUTOFF, FORCE_SYMM, CLUSTERS_PER_MUT_CUTOFF );

if ( print_mode )
    subplot(1,2,2)
end
show_2dmap( cluster_assignment_plot, structure, offset );

% for paper figure.
if ( print_mode );    
    figure(2)
    clf;
    subplot(1,1,1);
    show_2dmap( cluster_assignment_plot, structure, offset );
    outline_clusters( cluster_assignment_plot, cmap );
    cmap( 1,: ) = 1.0;
    colormap( cmap );
    xticklabel_rotate
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cluster_assignment  = get_cluster_filter( d, cluster_assignment, CLUSTER_CUTOFF, ZSCORE_SUM_CUTOFF, MAX_ZSCORE_CUTOFF, MEAN_ZSCORE_CUTOFF, NUM_MUT_CUTOFF, AVG_HITS_PER_MUT_CUTOFF, FORCE_SYMM, CLUSTERS_PER_MUT_CUTOFF );

NRES = length( d );

for i = 1:NRES
  n_clusters = length( unique( cluster_assignment( :, i ) ) );
  if ( n_clusters > CLUSTERS_PER_MUT_CUTOFF )
    cluster_assignment(:, i ) = 0;
  end
end

count = max( max( cluster_assignment) ) ;
for i = 1:count
  gp = find( cluster_assignment == i );
  if length(gp) > 0 & length(gp) < CLUSTER_CUTOFF
    cluster_assignment( gp ) = 0;
  end
end

for i = 1:count
  gp = find( cluster_assignment == i );
  if length(gp) > 0 & ( sum( d( gp ) ) < ZSCORE_SUM_CUTOFF | ...
			max( d( gp ) ) < MAX_ZSCORE_CUTOFF | ...
			mean( d(gp) ) < MEAN_ZSCORE_CUTOFF );
    cluster_assignment( gp ) = 0;
  end
end

[xgrid, ygrid ] = meshgrid( 1:NRES, 1:NRES );

for i = 1:count
  gp = find( cluster_assignment == i );
  if length(gp) > 0 
    num_mut_res = length( unique( xgrid( gp ) ) );
    if ( num_mut_res < NUM_MUT_CUTOFF )
      cluster_assignment( gp ) = 0;
    end
    if ( (length(gp) / num_mut_res ) >  AVG_HITS_PER_MUT_CUTOFF )
      cluster_assignment( gp ) = 0;
    end
    if FORCE_SYMM & ~check_for_symmetry( xgrid( gp ), ygrid( gp ) );
      cluster_assignment( gp ) = 0;
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = check_for_symmetry( x, y )
ok = 0;
if ~isempty( intersect( x, y ) )
  ok = 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cluster_assignment = check_neighbor(  cluster_assignment, i, j, ix, jx );

if ix > size( cluster_assignment, 1); return; end;
if jx > size( cluster_assignment, 2); return; end;
if ix < 1 ; return; end;
if jx < 1 ; return; end;

if cluster_assignment( ix, jx )  
  cluster_assignment = replace_cluster_assignment( cluster_assignment, ...
						   cluster_assignment(i,j), ...
						   cluster_assignment(ix, jx ) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   cluster_assignment = replace_cluster_assignment( cluster_assignment, a, b );

if (a == b); return; end;
gp = find( cluster_assignment == a );
cluster_assignment( gp ) = b;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outline_clusters( cluster_assignment_plot, cmap );
hold on

nres = size( cluster_assignment_plot, 1);
for i = 1:nres
  for j = 1:nres
    c = cluster_assignment_plot(i,j);
    if ( c > 0 )
      if ( i == 1 | ~cluster_assignment_plot(i-1,j)  )
	plot( [i-0.5 i-0.5], [j-0.5 j+0.5], 'color', 0.2 * cmap( c,: ) );
      end
      if ( i ==nres  | ~cluster_assignment_plot(i+1,j)  )
	plot( [i+0.5 i+0.5], [j-0.5 j+0.5], 'color', 0.2 * cmap( c,: ) );
      end
      if ( j == 1 | ~cluster_assignment_plot(i,j-1)  )
	plot( [i-0.5 i+0.5], [j-0.5 j-0.5], 'color', 0.2 * cmap( c,: ) );
      end
      if ( j == nres | ~cluster_assignment_plot(i,j+1)  )
	plot( [i-0.5 i+0.5], [j+0.5 j+0.5], 'color', 0.2 * cmap( c,: ) );
      end
    end
  end
end

hold off