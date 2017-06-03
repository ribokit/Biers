function show_recovered_helices( native_structure, structure, offset, sequence, Z );
% show_recovered_helices( native_structure, structure, Z );
%
%  Plotting function that shows helix recovery on both a 2D map (left)
%    and a secondary structure (right).
%
% INPUTS
%  native_structure = native (reference) structure in dot-parens notation
%  structure        = modeled structure in dot-parens notation
%  offset           = value to add to 1,2,... N to get conventional
%                        numbering
%  sequence         = sequence of RNA (use lowercase for positions to not display)
%  Z                = 2D map to display in panel 1
%

% figure out matches to native structure
native_stems  = parse_stems( native_structure );
modeled_stems = parse_stems( structure );
[native_recovered] = find_shared_stems( native_stems, modeled_stems );

% Colors for each stem
% cred, ... are attempt to reproduce "clarence" colors from M2-seq P4-P6
% depictions.
cred = [1 0.2 0.5]; ctan = [1,0.5,0.3]; cblue = [0.1 0.5 0.8]; cgold = [0.9, 0.5 0.3]; 
cgreen = [0.2 0.7 0.2]; cpurple = [1 0.2 1]; cpink = [1 0.8 0.8]; ccyan = [0.1 0.5 1];
conventional_colors = [cred;ctan;cblue;cgold;cgreen;cpurple; cpink; ccyan];
% fill out with MATLAB colors.
conventional_colors = [conventional_colors; jet(length(native_stems) - length(conventional_colors))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show Zscores in one pane
%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1);
if ~exist( 'Z' ) | isempty( Z )
    Z = zeros( length(structure),length(structure) );
    for i = 1:length( modeled_stems );
        stem = modeled_stems{i};
        for j = 1:size( stem, 1 );
            Z( stem(j,1), stem(j,2) ) = 1;
            Z( stem(j,2), stem(j,1) ) = 1;
        end
    end
    Z = 100*Z;
end
show_2dmap( Z, '', offset )
goodpos = find( upper( sequence ) == sequence );
axis( [min(goodpos), max(goodpos), min(goodpos), max(goodpos), ] );
set(gca,'fontsize',7);

% draw arrows to recovered domains on that 2D plot
hold on
stem_colors = conventional_colors( define_domains( native_stems ), : );
special_colors = {};
label_offset = length(structure)*(5/200);

for i = 1:length( native_stems )
    if ~native_recovered( i ); continue; end;
    stem = native_stems{i};
    x = mean(stem(:,1));
    y = mean(stem(:,2));
    % Could label helices by text if conventional_helix_names are supplied...
    %h1 = text( x+label_offset, y+label_offset,conventional_helix_names{i});
    %h2 = text( y+label_offset, x+label_offset,conventional_helix_names{i});
    %set([h1,h2], 'horizontalalign','center', 'fontweight','bold','fontsize',7,'color',stem_colors(i,:) );
    arrow( [x,y]+label_offset+1,[x,y]+2,'length',3,'width',1.5,'tipangle',30,'color',stem_colors(i,:) );
    arrow( [y,x]+label_offset+1,[y,x]+2,'length',3,'width',1.2,'tipangle',30,'color',stem_colors(i,:) );
    special_colors = [special_colors, stem_colors(i,:)];
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create varna fig highlighting recovered domains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2);
special_bps = native_stems( find( native_recovered ) );
for i = 1:length( special_bps ); special_bps{ i }  = special_bps{i} - min(goodpos ) + 1; end;
varna_fig([], sequence(goodpos), native_structure(goodpos), [], [], offset-min(goodpos)+1, special_bps, special_colors );
