function plot_and_save_sec_str_MOHCA(d, seq, str, file_name)

l = length(seq);

data_types = cell(1,l);
for i = 1:l
    data_types{i} = num2str(i);
end;
area_pred = generate_area_pred(seq,str,0,data_types,l);

figure();clf;
set_print_page(gcf,1);
image(d*35);colormap(1-gray()); axis equal; hold on;
for i = 1:l
    mark_points = find( area_pred(i, :)==1 );
    if ~isempty(mark_points);
        plot(mark_points + 1, i, 'o', 'MarkerEdgeColor', 'r');
    end;
    hold on;
end;
axis([1 l 1 l]);
title(file_name);

print(gcf, '-depsc2', '-loose', '-r300', file_name);
fprintf( ['Created: ', file_name, '\n'] );