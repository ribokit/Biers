function [ Z, mutpos, seqplot ] = output_Zscore_from_rdat( outfile, rdat_files, rdat_nomod, MEAN_EXPOSURE_CUTOFF, ZSCORE_OFFSET, APPLY_ZSCORE_OFFSET, ONLY_A_C, ignore_mut, print_stuff, mask_diag );
% [ Z, mutpos, seqplot ] = output_Zscore_from_rdat( outfile, rdat_files, rdat_nomod, MEAN_EXPOSURE_CUTOFF, ZSCORE_OFFSET, APPLY_ZSCORE_OFFSET, ONLY_A_C, ignore_mut, print_stuff, mask_diag );
%
% Z = [ reactivity value -  mean reactivity at residue across mutants]/ 
%                               standard-deviation at residue across mutants.
%
% Note that Z-scores are multiplied by -1.0 for output.
%
%  outfile    = name of an output text file
%  rdat_files = name of mutate-and-map file in RDAT format. [or several files if adding values]
%
% Optional:
%  rdat_nomod = name of mutate-and-map file, control rxn without modiifer. 
%                Used for filtering. Give [] if you don't have one. (Default: [])
%  MEAN_EXPOSURE_CUTOFF = cutoff on mean peak intensity -- residues with mean reactivity
%                         above this cutoff will be considered already exposed, and not used.
%                         note that this cutoff applies to values after 
%                         normalizing within each mutant (Default: 1.0)
%  ZSCORE_OFFSET = offset Z-score by this value, and only output values > 0.0. (Default: 0.0)
%  APPLY_ZSCORE_OFFSET = apply the offset (Default: 1)
%  ONLY_A_C = only calculate Z-scores at A or C residues (may be useful for DMS data).
%               (Default: 0)
%  ignore_mut  = ignore variants with mutations at the specified positions. (Default: [])
%  print_stuff = verbose output (Default: 0)
%  mask_diag   = remove data +/- this value from diag before Z-score compute (Default: -1 means no masking)
%
% Originally part of HiTRACE tools, now in Biers with other RNA secondary structure modeling commands.
% 
% (C) R. Das, Stanford University, 2010-2013, 2017
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'MEAN_EXPOSURE_CUTOFF','var' ); MEAN_EXPOSURE_CUTOFF = 1.0; end;
if ~exist( 'APPLY_ZSCORE_OFFSET','var' ); APPLY_ZSCORE_OFFSET = 1; end;
if ~exist( 'ZSCORE_OFFSET','var' ); ZSCORE_OFFSET = 0.0; end;
if ~exist( 'ONLY_A_C','var' ); ONLY_A_C = 0; end;
if ~exist( 'print_stuff','var' ); print_stuff = 0; end
if ~exist( 'rdat_nomod','var' ) || isempty( rdat_nomod ); rdat_nomod = []; end
if ~exist( 'ignore_mut','var' ); ignore_mut = []; end

Z = []; mutpos = [];seqplot = [];
if isempty( rdat_files ); return; end;

if ~iscell( rdat_files ); rdat_files = { rdat_files }; end;

NRES = 0;
if ~isempty( rdat_nomod); 
  if ischar( rdat_nomod )
    rdat_nomod = read_rdat_file( rdat_nomod);
  end
  NRES = length( rdat_nomod.sequence ); 
  Z_sum = zeros( NRES, NRES );
end;
mut_weights_sum = zeros( 1, NRES );

for i = 1:length( rdat_files )

  if ( ischar( rdat_files{i} ) )
      d = read_rdat_file( rdat_files{i} );
  else
      d = rdat_files{i};
      assert( isobject( d ) );
  end
  
  if ~exist( 'mask_diag','var' );
      % look for a stripe of 1.00's -- if
      if ( length( find( abs(d.reactivity - 1.00) < 1e-5 ) ) > 10 )
          mask_diag = 5;
          fprintf( 'Looks like an M2-seq data set based on ones on the diagonal. Masking positions within %d of diagonal.\n', mask_diag )
      else
          mask_diag = -1;
      end
  end

  if ( NRES == 0 )
    NRES = length(d.sequence);
    Z_sum = zeros( NRES, NRES );
    mut_weights_sum = zeros( 1, NRES );
  elseif ( length(d.sequence) ~= NRES ) 
    fprintf( 'WARNING! WARNING! NRES problem! %s\n', rdat_files{i} );
  end


  [ Z, mutpos ] = get_Zscore_and_apply_filter( d, rdat_nomod, MEAN_EXPOSURE_CUTOFF, ZSCORE_OFFSET, APPLY_ZSCORE_OFFSET, ONLY_A_C, ignore_mut, mask_diag );
  % just a check
  if ~isempty(strfind( d.name, 'P4P6'));    Z( 176-d.offset, : ) = 0.0;   end;

  Z_sum = Z_sum + Z; % + smooth2d( Z );  

  mut_weights = ( sum( Z, 1 ) < 0 );
  mut_weights_sum = mut_weights_sum + mut_weights;

end


for i = 1:NRES
  if ( mut_weights_sum(i) > 0 )
    Z(:,i) = Z_sum(:,i)/mut_weights_sum(i);    
  end
end

seqplot = [1:length( d.sequence )] + d.offset;
plot_and_save(Z, seqplot,  outfile, print_stuff );
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Zscore_full, mutpos ] = get_Zscore_and_apply_filter( d, rdat_nomod, MEAN_EXPOSURE_CUTOFF, ZSCORE_OFFSET, APPLY_ZSCORE_OFFSET, ONLY_A_C, ignore_mut, mask_diag );

Zscore = [];

a_input = d.reactivity;

a = a_input;
a = mask_diagonal( a, mask_diag );

normbins = [ 10: size( a,1)-10 ];
a = quick_norm( a, normbins );

for i = 1:size( a, 1 )
  Zscore(i,:) =  ( a(i,:) - mean( a(i,:)) )/ std( a(i,:),0,2);
  mean_exposure(i) = mean( a(i,:) );
end

% special case -- what if some columns are different experiments.
diff_ex_cols = [];
for j = 1:size( a, 2 )
  data_annotation = d.data_annotations{j};
  for k = 1:length( data_annotation )
    if ( strfind( data_annotation{k}, 'differentExperiment' ) )
      diff_ex_cols = [ diff_ex_cols, j ];
      break;
    end
  end
end
%if ~isempty( 'diff_ex_cols' )
%  diff_ex_cols
%end

for i = 1:size( a, 1 )
  Zscore(i,diff_ex_cols) =  ( a(i,diff_ex_cols) - mean( a(i,diff_ex_cols)) )/ std( a(i,diff_ex_cols),0,2);
end

% Don't include rows that have uniformly high chemical modification --
% mutate/map looks for 'release' of protected residues.
for i = 1:size( a, 1 )
  if ( mean_exposure(i) > MEAN_EXPOSURE_CUTOFF  )
    Zscore(i,:) = 0.0;
  end
end

if ( ONLY_A_C )
  for i = 1:size( a, 1 )
    seqchar =   d.sequence( d.seqpos( i ) - d.offset );
    if ( seqchar ~='A' & seqchar ~='C' )
      Zscore(i,:) = 0.0;
    end  
  end
end

% Filter to make sure Z-score is above some cutoff.
if APPLY_ZSCORE_OFFSET
  Zscore = max( ( Zscore - ZSCORE_OFFSET ), 0.0 );
end

% filter for outlier columns (typically bad RNA in that well...)
%MEAN_Z_FILTER = 1;
%if MEAN_Z_FILTER
%  mean_Z = mean( abs( Zscore ) );
%
%  for j = 1:size( a_input, 2 );
%    if ( mean_Z(:,j) > mean(mean_Z) * 4 )
%      Zscore(:,j) = 0.0;
%    end
%  end
%end

for j = 1:size( a, 2 )
  data_annotation = d.data_annotations{j};
  for k = 1:length( data_annotation )
    if ( contains( data_annotation{k}, 'badQuality' ) )
      Zscore(:,j) = 0.0;
    end
  end
end


% If 'nomod' (i.e., background) measurement is given, use it to get rid of strong background spots.
Z_NOMOD_CUTOFF = 6.0;
if ~isempty( rdat_nomod )  

  b = rdat_nomod.reactivity;
  if ( size( b, 1 ) ~= size( a, 1 ) )
    fprintf( 'Hey, rdat_nomod needs to have the same dimensions as d !!!\n' ); return;
  end
  b = quick_norm( b, normbins );

  for i = 1:size( b, 1 )
    Zscore_nomod(i,:) =  ( b(i,:) - mean( b(i,:)) )/ std( b(i,:),0,2);
  end  

  for i = 1:size(b,1)
    for j = 1:size(b,2)
      if Zscore_nomod(i,j) > Z_NOMOD_CUTOFF
	Zscore( i,j ) = 0.0;
      end
    end
  end

end

% Need to reshape into NRES by NRES matrix.
% Convert factor
ZSCORE_SCALING = -1.0;

if ~exist('d.mutpos','var') || length( d.mutpos ) == 0
  mutpos = generate_mutpos_from_rdat(d.data_annotations,d);
end

NRES = length(d.sequence);
Zscore_full = zeros( NRES, NRES );
pos1 = d.seqpos - d.offset;
pos2 = mutpos; % - d.offset;

% look for first library, as marked by MutPos. Note that first column is assumed to be wt.
for i = 2 : length( mutpos )
  if (i < length( mutpos )) & ( isnan( mutpos(i+1) ) || mutpos(i) > mutpos(i+1) ); 
    break; 
  end
end
gp = [2:i];

% Remove mutants that are inputted in "ignore_mut"
for i = 1:length( ignore_mut )
  gp = setdiff( gp,   find( mutpos == ignore_mut(i)) );
end
Zscore_full( pos1, pos2(gp) ) = Zscore(:,gp) * ZSCORE_SCALING;

mutpos = mutpos( gp );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_and_save(Z, seqplot,  outfile, print_stuff );

if ~exist( 'print_stuff' ); print_stuff = 0; end
%image( seqplot, seqplot,  (-Z' - 1.5 ) * 64 );
image( seqplot, seqplot,  (-Z') * 64 );

hold on
plot( seqplot, seqplot, 'k' ); 
hold off
colormap( 1 - gray(100) );
h=title( outfile );set(h,'interpreter','none');
xlabel( 'seqpos');
ylabel( 'mutpos' );


if length( outfile ) > 0
    save( outfile, '-ascii','Z');
    if ( print_stuff )
    set(gcf,'PaperPositionMode','Auto');
    %set(gcf,'Position',[0, 0, 600, 600]);
        eps_file = [outfile, '.eps' ];
        fprintf( 'Outputting to postscript file: %s\n', eps_file );
        print( eps_file, '-depsc2', '-tiff' );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  d = fill_mutpos( d )

% need to figure out where mutations are based on tags like "mutation:G64C" in data_annotation
d.mutpos = [];
for k = 1:length( d.data_annotations )
  d.mutpos(k) = NaN;
  data_annotation = d.data_annotations{k};
  for m = 1:length( data_annotation )
    c = str2cell( data_annotation{m},':' );
    if length(c)> 0 & strcmp( c{1}, 'mutation' )
      num = str2num( remove_AGCTU( c{2} ) );
      if length(num)>0;  d.mutpos(k) =  num; end;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = str2cell(s, delim)

if ~exist( 'delim' ) 
  tabchar = sprintf( '\t' );
  if ~isempty(strfind( s, tabchar ) ) 
    delim = tabchar;
  else
    delim = ' '; 
  end
end;

rest = s;
i = 1;
c = {};
while length(rest)
    [t, rest] = strtok(rest, delim);
    c{i} = t;
    i = i + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = remove_AGCTU( x )
x = strrep( x, 'A', '');
x = strrep( x, 'G', '');
x = strrep( x, 'C', '');
x = strrep( x, 'T', '');
x = strrep( x, 'U', '');

