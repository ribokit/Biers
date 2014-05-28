function mohcaplot_biers( D, seqpos, ligpos, titl, ticksize, save_path, secstr, c )

% Plots 2D maps in MOHCA style
%
% INPUTS:
%       D         req  = matrix of data to be plotted
%       seqpos    opt  = x-axis values, RT stop positions (enter '' for default, 1 to length of x-axis data in D) 
%       ligpos    opt  = y-axis values, ligation positions (enter '' for default, 1 to length of x-axis data in D)
%       titl      opt  = desired plot title, string (enter '' for default, no title) 
%       ticksize  opt  = font size of tick labels (default 25, enter '' for default) 
%       save_path opt  = path to save file (including filename) (if none, enter '') 
%       secstr    opt  = cell array with {sequence, secstr, offset, data_types, numlanes} 
%       c         opt  = colorcode (1 for jet, 2 for gray)
%
% Clarence Cheng, 2014
%

if ~exist( 'seqpos', 'var' ) || isempty( seqpos ); seqpos = 1:size(D,2); end
if ~exist( 'ligpos', 'var' ) || isempty( ligpos ); ligpos = 1:size(D,1); end
if ~exist( 'ticksize', 'var' ) || isempty( ticksize ); ticksize = 12; end
if ~exist( 'titl', 'var' ); titl = ''; end
if ~exist( 'c', 'var') || isempty( c ); c = 1; end;
% if ~exist( 'contours', 'var' ); contours = ''; end

% Make plot
figure;
set(gcf, 'PaperPositionMode', 'Manual','PaperOrientation', 'Landscape','PaperPosition', [-0.65 0.15 12 8],'color','white','Position', [0 0 800 600]);
image( seqpos, ligpos, 50 * D' ); hold on;
axis image;
if c == 1; colormap( jet); else colormap( 1-gray ); end;

% Label x and y axes
gp = find( mod(seqpos,10) == 0 );
set(gca,'xtick',seqpos(gp) )
gp = find( mod(ligpos,10) == 0 );
set(gca,'ytick',ligpos(gp) )
set(gca,'TickDir','out');
set(gca,'xgrid','on','ygrid','on','fonts',ticksize);
xlabel( 'Reverse transcription stop position [5'']','fontsize',15,'fontweight','bold' );
ylabel( 'Cleaved and ligated position [3'']','fontsize',15,'fontweight','bold' );
hold on;

% Add title
title( titl, 'fonts', 15, 'fontw', 'bold' );

% Rotate xticklabels and reposition
xticklabel = get(gca,'XTickLabel');
set(gca,'XTickLabel','');
hxLabel=get(gca,'XLabel');
set(hxLabel,'Units','data');
xLabelPosition=get(hxLabel,'Position');
y=xLabelPosition(2) - 7;
XTick=str2num(xticklabel)+1;
y=repmat(y,length(XTick),1);
fs = get(gca,'fonts');
hText=text(XTick,y,xticklabel,'fonts',ticksize);
set(hText,'Rotation',90,'HorizontalAlignment','right');
xlab=get(gca,'XLabel');
set(xlab,'Position',get(xlab,'Position') + [0 7 0]);

% Make colorbar legend
hc = colorbar('location','eastoutside');
hcm = max(get(hc,'YLim'));
set(hc,'YTick',[0.5 hcm-0.5]);
set(hc,'YTickLabel',{'0.0','1.0'});
hcp = get(hc,'pos');
pos = get(gca,'pos');
%set(hc,'position',[hcp(1)*0.92 hcp(2) hcp(3)*0.5 pos(4)*0.25],'fonts',15,'fontw','bold');

% Overlay secondary structure
if exist( 'secstr', 'var' )
    seq = secstr{1};
    str = secstr{2};
    data_types = cell(1,length(secstr{1}));
    for i = 1:length(secstr{1})
        data_types{i} = num2str(i);
    end
    area_pred = triu(generate_area_pred(seq, str, 0, data_types, length(secstr{1})));
        % in future, use sequence and secstruct from rdat and crop to correct size
    [x_pred, y_pred] = find(area_pred);
    x_pred = x_pred - 1 + seqpos(1);
    y_pred = y_pred - 1 + ligpos(1);
    plot(x_pred, y_pred, '.', 'color', [1 0 1]);
end

% Save figure
if exist( 'save_path', 'var' )
    if ~isempty( save_path )
        print( gcf, '-depsc2', '-loose', '-r300', save_path);
        fprintf( ['Created: ', save_path, '\n'] );
    end
end