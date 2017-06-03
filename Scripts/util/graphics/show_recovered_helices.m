function show_recovered_helices( native_structure, structure, Z );

% figure out matches to native structure
native_stems = parse_stems( r.structure );
modeled_stems = parse_stems( structure );

[native_recovered] = find_shared_stems( native_stems, modeled_stems );

subplot(1,2,1);
show_2dmap( -20*min(Zscores{i},0), '', r.offset )
goodpos = find( upper( r.sequence ) == r.sequence );
axis( [min(goodpos), max(goodpos), min(goodpos), max(goodpos), ] );
set(gca,'fontsize',7);

% draw arrows to recovered domains

hold on
stem_colors = conventional_colors( define_domains( native_stems ), : );
special_colors = {};
label_offset = length(r.sequence)*(5/200);
for i = 1:length( native_stems )
    if ~native_recovered( i ); continue; end;
    stem = native_stems{i};
    x = mean(stem(:,1));
    y = mean(stem(:,2));
    % do this if conventional_helix_names are supplied...
    %h1 = text( x+label_offset, y+label_offset,conventional_helix_names{i});
    %h2 = text( y+label_offset, x+label_offset,conventional_helix_names{i});
    %set([h1,h2], 'horizontalalign','center', 'fontweight','bold','fontsize',7,'color',stem_colors(i,:) );
    arrow( [x,y]+label_offset+1,[x,y]+2,'length',3,'width',1.5,'tipangle',30,'color',stem_colors(i,:) );
    arrow( [y,x]+label_offset+1,[y,x]+2,'length',3,'width',1.2,'tipangle',30,'color',stem_colors(i,:) );
    special_colors = [special_colors, stem_colors(i,:)];
end
hold off

subplot(1,2,2);
% create varna fig highlighting recovered domains.
special_bps = native_stems( find( native_recovered ) );
for i = 1:length( special_bps ); special_bps{ i }  = special_bps{i} - min(goodpos ) + 1; end;
varna_fig([], r.sequence(goodpos), r.structure(goodpos), [], [], r.offset-min(goodpos)+1, special_bps, special_colors );
pause;