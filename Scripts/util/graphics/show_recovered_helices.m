function show_recovered_helices( native_structure, structure, offset, sequence, Z, nt_labels, which_subplots, tag );
% show_recovered_helices( native_structure, structure, offset, sequence, Z, nt_labels, which_subplots, tag );
%
%  Plotting function that shows helix recovery on both a 2D map (left)
%    and a secondary structure (right).
%
%  Tries to be smart about where to apply nucleotide labels (nt_labels).
%
%  NOTE: Smooths 2d Z-scores for visualization! 
%
% INPUTS
%  native_structure = native (reference) structure in dot-parens notation
%  structure        = modeled structure in dot-parens notation
%  offset           = value to add to 1,2,... N to get conventional
%                        numbering
%  sequence         = sequence of RNA (use lowercase for positions to not display)
%  Z                = 2D map to display in panel 1
%  nt_labels        = cell of labels for particular nts in helices, like { {149,'P5b'},... };
%  which_subplots   = two triplets defining MATLAB subplot locations (default: [] is left and right panel)
%
% (C) R. Das, Stanford University, June 2017

% figure out matches to native structure
native_stems  = parse_stems( native_structure );
modeled_stems = parse_stems( structure );
[native_recovered] = find_shared_stems( native_stems, modeled_stems );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create 2D plot with annotations for detected helices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist( 'which_subplots', 'var' ) | isempty( which_subplots ); which_subplots = { [1,2,1], [1,2,2] }; end;
subplot(which_subplots{1}(1),which_subplots{1}(2),which_subplots{1}(3));
if ~exist( 'Z' ) | isempty( Z )
    Z = zeros( length(structure),length(structure) );
    for i = 1:length( modeled_stems );
        stem = modeled_stems{i};
        for j = 1:size( stem, 1 );
            Z( stem(j,1), stem(j,2) ) = 1;
            Z( stem(j,2), stem(j,1) ) = 1;
        end
    end
else
    Z = smooth2d( Z );
end
show_2dmap( Z+0.5, '', offset, -2 )
set(gca,'xgrid','off','ygrid','off' );
goodpos = find( upper( sequence ) == sequence );
axis( [min(goodpos), max(goodpos), min(goodpos), max(goodpos), ] );
set(gca,'fontsize',7);


hold on
if ~exist( 'nt_labels' ) | isempty(nt_labels); nt_labels = {}; end;
stem_colors = get_conventional_stem_colors( native_stems );
stem_labels = assign_stem_labels( native_stems, nt_labels, offset );
special_colors = {};
label_offset = length(structure)*(5/200);
stem_labels_to_display = {};
for i = 1:length( native_stems )
    if ~native_recovered( i ); continue; end;
    stem = native_stems{i};

    % Label helices by text if conventional_helix_names are supplied... 
    % figure out positions based on Z-scores -- and accumulate across
    % any stems that have the same label.
    stem_label = stem_labels{i};
    if length( stem_label ) > 0
        idx = find( strcmp( stem_labels_to_display, stem_label ) );
        if ( isempty(idx) );
            idx = length( stem_labels_to_display ) + 1;
            stem_labels_to_display{idx} = stem_labels{i};
            stem_label_colors(idx,:) = stem_colors(i,:);
            wtT( idx ) = 0;        pos_wtdT(:,idx) = [0,0];
            wtB( idx ) = 0;        pos_wtdB(:,idx) = [0,0];
        end;
        for n = 1:size( stem, 1 )
            wtT(idx ) = wtT(idx) + Z( stem(n,1),stem(n,2) );
            pos_wtdT(:,idx) = pos_wtdT(:,idx) + Z( stem(n,1),stem(n,2) )* stem(n,:)';
            wtB(idx ) = wtB(idx) + Z( stem(n,2),stem(n,1) );
            pos_wtdB(:,idx) = pos_wtdB(:,idx) + Z( stem(n,2),stem(n,1) )* fliplr(stem(n,:)');
        end
    else
        x = mean(stem(:,1));
        y = mean(stem(:,2));
        arrow( [x,y]+label_offset+1,[x,y]+2,'length',3,'width',1.5,'tipangle',30,'color',stem_colors(i,:) );
        arrow( [y,x]+label_offset+1,[y,x]+2,'length',3,'width',1.2,'tipangle',30,'color',stem_colors(i,:) );
    end
    special_colors = [special_colors, stem_colors(i,:)];
end

extra_annotations = {};
for i = 1:length( stem_labels_to_display )
    x = pos_wtdT(1,i)/wtT(i);    y = pos_wtdT(2,i)/wtT(i);
    h1 = text( x+label_offset, y+label_offset,stem_labels_to_display{i});
    x = pos_wtdB(1,i)/wtB(i);    y = pos_wtdB(2,i)/wtB(i);
    h2 = text( y+label_offset, x+label_offset,stem_labels_to_display{i});
    set([h1,h2], 'horizontalalign','center', 'fontweight','bold','fontsize',9,'color',stem_label_colors(i,:) );
    extra_annotations{i} = {stem_labels_to_display{i}, round( x ), stem_label_colors(i,:) };
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create varna fig highlighting recovered domains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(which_subplots{2}(1),which_subplots{2}(2),which_subplots{2}(3));
special_bps = native_stems( find( native_recovered ) );
for i = 1:length( special_bps ); special_bps{ i }  = special_bps{i} - min(goodpos ) + 1; end;
for i = 1:length( extra_annotations ); extra_annotations{i}{2}  = extra_annotations{i}{2} - min(goodpos ) + 1; end;

% PNG (will show up in MATLAB);
varna_fig([], sequence(goodpos), native_structure(goodpos), [], [], offset+min(goodpos)-1, special_bps, special_colors, [], [], extra_annotations );
if ~exist( 'tag','var' ) tag = 'show_recovered_helices'; end;
% EPS -- editable in illustrator
varna_fig([tag,'.eps'], sequence(goodpos), native_structure(goodpos), [], [], offset+min(goodpos)-1, special_bps, special_colors, [], [], extra_annotations );
%varna_fig([tag,'.svg'], sequence(goodpos), native_structure(goodpos), [], [], offset-min(goodpos)+1, special_bps, special_colors, [], [], extra_annotations );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stem_labels, domain_labels] = assign_stem_labels( native_stems, nt_labels, offset );
stem_labels = {}; domain_labels = {};

for i = 1:length( native_stems ); stem_labels{i} = ''; end;

 % domain labels 'inherit' from stem labels.
domains = define_domains( native_stems );
for i = 1:max( domains ); domain_labels{i} = ''; end;

% figure out which stem each label sits in.
% that will also be the 'default' label for the stem's domain.
for j = 1:length( nt_labels )
    nt = nt_labels{j}{1};
    label = nt_labels{j}{2};
    for i = 1:length( native_stems )
        stem = native_stems{i};
        if ~isempty( find( (stem+offset) == nt ) )
            stem_labels{i} = label;
            domain_labels{ domains(i) } = label;
            break;
        end
    end
end

% give labels to other stems:
for i = 1:length( native_stems )
    stem = native_stems{i};
    if ( length( stem_labels{i} ) == 0 )
        %assert( length( domain_labels{ domains(i) } ) > 0 );
        stem_labels{i} = domain_labels{ domains(i) } ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x_smooth = smooth2d( x, niter );
% x_smooth = smooth2d( x, niter );
if ~exist( 'niter' )  niter = 2; end
B = [ 0 0.1 0; 0.1 0.6 0.1; 0 0.1 0];
x_smooth = conv2( x, B, 'same' );

for  n = 1:niter-1;
    x_smooth = conv2( x_smooth, B, 'same' );
end

