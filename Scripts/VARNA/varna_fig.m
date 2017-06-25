function varna_fig(filename, sequence, structure, DATA, colorscheme, offset, special_base_pairs, special_colors, bpp_values, bpp_anchor_bases, extra_annotations )
% VARNA_FIG: Create image file with secondary structure -- double click to get VARNA visualizer in a web browser
%
%  varna_fig(filename,sequence,structure,DATA,colorscheme,offset,special_base_pairs,special_colors, bpp_values, bpp_anchor_bases)
%
% filename  = output filename [end in '.eps','.png','.svg']. PNG files will
%              display in MATLAB. Default is [], which is interpreted as
%              /tmp/tmp.png.
% sequence  = RNA sequence
% structure = structure in dot/bracket notation. length should match sequence.
%
% optional input parameters:
% DATA               = data for coloring residues. will only use values between 0 and 1,
%                        so normalize ahead of time. Give [] if no data.
% colorscheme        = blue/white/red (0), another blue/white/red (1), white/orange/red (2)
% offset             = number to add to 1,2,...N to get conventional numbering
% special_base_pairs = extra pairs of residues that should have connecting lines.
% special_colors     = colors of connecting lines for special_base_pairs.
% bpp_values         = numbers that will show up on helices as percentages
% bpp_anchor_bases   = where those bpp_values will show up.
% extra_annotations  = cell of label-resnum-rgb triplets, like 
%                            { {'P5b',150,[0.5 0.5 0.5]}, ... }.
%
% (C) R. Das 2011, 2017
% (C) C.C. VanLang, P. Cordero 2010
% (C) S. Tian, 2016

if nargin < 3; help(mfilename); return; end;
if isempty( filename ); filename = '/tmp/tmp.png'; end;
if ~exist('colorscheme', 'var'); colorscheme = 1; end;

if ~isempty(DATA)
    % DATA Prep
    %%%reactivity=hSHAPE(DATA);
    reactivity = DATA;
    
    reactivity = max(reactivity, 0);
    reactivity = min(reactivity, 2.0);
    %reactivity=100*reactivity;
    
    graypoints = find(DATA == -999 | isnan(DATA));
    reactivity(graypoints) = -0.01;
end;

is_tmp_file = ~isempty(find(strfind( filename, '/tmp/')==1));

if exist('special_base_pairs', 'var');
    if length(special_base_pairs) ~= length(special_colors);
        fprintf('Must specify a special_color for each special_base_pair set\n');
    end;
end

if length( filename ) < 5 | ~strcmp( filename( end-4:end), '.html') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now standard mode -- use command-line interface into VARNA.jar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    [VARNA_DIR, VARNA_JAR] = get_varna();
    if ~exist( [VARNA_DIR,VARNA_JAR], 'file' ); 
        fprintf( ['Cannot find: ',VARNA_DIR,VARNA_JAR,'\n' ] );
        scripts_dir = fileparts( fileparts( which( 'varna_fig.m') ) );
        fprintf( ['Please update ',scripts_dir,'/get_varna.m based on ',scripts_dir,'/get_varna.m.example\n' ]);
        return;
    end;
    
    command = ['java -cp ',VARNA_DIR,VARNA_JAR,' fr.orsay.lri.varna.applications.VARNAcmd', ' \\\n' ];
    command = [command, ' -sequenceDBN ', sequence, ' \\\n'];
    command = [command, ' -structureDBN "', structure,'"', ' \\\n'];
    command = [command, ' -algorithm radiate', ' \\\n'];
    
    if ~isempty(DATA);
        command = [command, ' -colorMap "'];
        for i = 1:length(reactivity);
            command = [command,sprintf('%6.3f;', reactivity(i)) ];
        end;
        command = [ command, '"', ' \\\n'];
    end
    
     if exist('reactivity', 'var');
         switch colorscheme
             case 0 % previous default
                 command = [ command, ' -colorMapStyle "0:#0000FF;10:#0000FF;40:#FFFFFF;60:#FFFFFF;90:#FF0000;100:#FF0000"', ' \\\n' ];
             case 1 % new default
                 command = [ command, ' -colorMapStyle "-0.01:#B0B0B0;0:#0000FF;1:#FFFFFF;2:#FF0000"', ' \\\n' ];
             case 2 % white orange to red
                 % slight pain because of VARNA rescaling:
                 if (sum( reactivity < 0 ) > 0)
                     command = [ command, ' -colorMapStyle "-0.001:#C0C0C0,0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000"', ' \\\n' ];
                 else
                     command = [ command, ' -colorMapStyle "0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000"', ' \\\n' ];
                 end
            case 3 % blue to white to yellow
                % slight pain because of VARNA rescaling:
                if (sum( reactivity < 0 ) > 0);
                    command = [ command, ' -colorMapStyle "-0.001:#B0B0B0,0:#0000FF;0.5:#FFFFFF;1:#FFFF00"', ' \\\n' ];
                else
                    command = [ command, ' -colorMapStyle "0:#FFFFFF;,0.5:#2222FF;1:#0000FF"', ' \\\n' ];
                end;
        end;
    end;
    
    command = [command, ' -bpStyle lw', ' \\\n'];
    command = [command, ' -baseInner "#FFFFFF"', ' \\\n'];
    command = [command, ' -baseOutline "#FFFFFF"', ' \\\n'];
    command = [command, ' -bp "#000000"', ' \\\n'];
    command = [command, ' -spaceBetweenBases 0.6', ' \\\n'];
    command = [command, ' -flat false', ' \\\n'];
    if ~is_tmp_file
        titlename = filename;
        dots = strfind(titlename,'.');
        if ~isempty(dots);  titlename = titlename( 1: dots(end)-1 ); end;
        command = [command, ' -title ',titlename, ' \\\n' ];
        command = [command, ' -titleColor "#000000"', ' \\\n'];
        command = [command, ' -titleSize 20', ' \\\n'];
    end
    command = [command, ' -colorMapCaption Reactivity', ' \\\n'];
        
    if exist('special_base_pairs', 'var') & length( special_base_pairs ) > 0 & length( special_base_pairs{1} ) > 0;
        command = [command, ' -auxBPs "'];         
        for q = 1:length(special_base_pairs)
            special_base_pair_set = special_base_pairs{q};
            hex_color = convert_rgb_to_hexadecimal(special_colors{q});
            for k = 1:size(special_base_pair_set, 1);
                command = [command, sprintf('(%d,%d):thickness=3,color=#%6s;', special_base_pair_set(k, 1), special_base_pair_set(k, 2), hex_color) ];
            end;
        end;    
        command = [command, '"', ' \\\n'];         
    end
    
    if exist('offset', 'var');
        command = [command, ' -baseNum "#FFFFFF"', ' \\\n'];
        command = [command, ' -periodNum 1000', ' \\\n'];
    end;

    if (exist('offset', 'var') | exist('bpp_values', 'var')) | exist( 'extra_annotations', 'var' )
         
        command = [command, ' -annotations "'];
        if exist('offset', 'var');
            PERIOD = 20;
            for i = 1:length(sequence)
                if (mod(i+offset, PERIOD) == 0)
                    command = [command, sprintf( '%d:type=B,anchor=%d,color=#000000,size=8;', i+offset, i )] ;
                end;
            end;
        end;

        if exist('bpp_values', 'var') && ~isempty(bpp_values);
            if length(bpp_values) ~= length(bpp_anchor_bases);
                fprintf('Must specify a bpp_anchor_base for each bpp_value \n');
            end;
            for i = 1:length(bpp_values)
                command = [command, sprintf('%3.0f%%:type=L,anchor=%d,color=#009000,size=9;', 100 * bpp_values(i), bpp_anchor_bases(i) )];
            end;
        end;        

        if exist('extra_annotations', 'var') && ~isempty(extra_annotations);
            for i = 1:length(extra_annotations)
                command = [command, sprintf('%s:type=B,anchor=%d,color=#%6s,size=12;', extra_annotations{i}{1}, ...
                    extra_annotations{i}{2},convert_rgb_to_hexadecimal(extra_annotations{i}{3}) )];
            end;
        end;        
        
        command = [command,'"', ' \\\n'];
    end

    command = [command, ' -flat true', ' \\\n' ];
    command = [command, ' -resolution 4.0', ' \\\n'];
    command_without_output = command;
    command = [command, ' -o ',filename];

    fprintf( [strrep(command, '%:', '%%:'), '\n'] );
    returncode = system( strrep(command, ' \\\n', '') );

    if ( returncode == 0 & ~is_tmp_file & system( 'which open > /dev/null' ) == 0 ) ;
        system( ['open ', filename ] ); 
    end;
    
    % show .png in matlab
    [dirname,basename,ext] = fileparts( filename );
    if ( strcmp(ext,'png') )  
        fprintf( 'Displaying PNG in MATLAB\n' );
        imshow( imread( filename ) );
    end

    % will display interactive VARNA
    if ( length( ext ) == 0 ) returncode = system( [strrep(command_without_output, ' \\\n', ''), ' &' ] ); end;

else
    
    % Old HTML-based output -- as of 2017, no Web browsers allow for Java to run in browser, so this should be deprecated soon.
    fid = fopen(filename,'w');
    
    vers = 1;
    [VARNA_DIR, VARNA_JAR] = get_varna();
    
    fprintf(fid, '%s\n', '<HTML><HEAD><META http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></HEAD><BODY>');
    fprintf(fid, '%s%s%s%s%s\n', '<APPLET code="VARNA.class" codebase="', VARNA_DIR, '" archive="', VARNA_JAR, '" width="1200" height="1200">');
    fprintf(fid, '%s%s%s\n', '<PARAM name="sequenceDBN"  value="', sequence, '"/>');
    fprintf(fid, '%s%s%s\n', '<PARAM name="structureDBN" value="', structure, '"/>');
    % fprintf(fid, '%s\n', '<PARAM name="algorithm" value="naview" />');
    fprintf(fid, '%s\n', '<PARAM name="algorithm" value="radiate" />');
    
    if ~isempty(DATA);
        fprintf(fid, '%s', '<param name="colorMap" value="');
        for i = 1:length(reactivity);
            fprintf(fid,' %6.3f%s', reactivity(i), ',');
        end;
        fprintf(fid, '%s\n','"/>');
    end;
    
    if exist('reactivity', 'var');
        switch colorscheme
            case 0 % previous default
                fprintf(fid, '%s\n', '<param name="colorMapStyle" value="0:#0000FF;10:#0000FF;40:#FFFFFF;60:#FFFFFF;90:#FF0000;100:#FF0000" />');
            case 1 % new default
                fprintf(fid, '%s\n', '<param name="colorMapStyle" value="-0.01:#B0B0B0;0:#0000FF;1:#FFFFFF;2:#FF0000" />');
            case 2 % white orange to red
                % slight pain because of VARNA rescaling:
                if (sum( reactivity < 0 ) > 0)
                    fprintf(fid, '%s\n', '<param name="colorMapStyle" value="-0.001:#C0C0C0,0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000" />');
                else
                    fprintf(fid, '%s\n', '<param name="colorMapStyle" value="0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000" />');
                end
            case 3 % blue to white to yellow
                % slight pain because of VARNA rescaling:
                if (sum( reactivity < 0 ) > 0);
                    fprintf(fid, '%s\n', '<param name="colorMapStyle" value="-0.001:#B0B0B0,0:#0000FF;0.5:#FFFFFF;1:#FFFF00" />');
                else
                    fprintf(fid, '%s\n', '<param name="colorMapStyle" value="0:#FFFFFF;,0.5:#2222FF;1:#0000FF" />');
                end;
        end;
    end;
    
    fprintf(fid, '<param name="bpStyle" value="lw" />\n');
    fprintf(fid, '<param name="baseInner" value="#FFFFFF" />\n');
    fprintf(fid, '<param name="baseOutline" value="#FFFFFF" />\n');
    fprintf(fid, '<param name="bp" value="#000000" />\n');
    fprintf(fid, '<param name="spaceBetweenBases" value="0.6" />\n');
    fprintf(fid, '<param name="flat" value="false" />\n');
    fprintf(fid, '%s%s%s\n', '<param name="title" value="', strrep(filename, '.html', ''), '" />\n');
    fprintf(fid, '<param name="titleColor" value="#000000" />\n');
    fprintf(fid, '<param name="titleSize" value="20" />\n');
    fprintf(fid, '<param name="colorMapCaption" value="Reactivity" />\n');
    
    % fprintf(fid, '<param name="colorMapMax" value="2" />\n');
    % fprintf(fid, '<param name="colorMapMin" value="0" />\n');
    
    if exist('special_base_pairs', 'var');
        if length(special_base_pairs) ~= length(special_colors);
            fprintf('Must specify a special_color for each special_base_pair set\n');
        end;
        
        fprintf(fid, '<param name="auxBPs" value="');
        
        for q = 1:length(special_base_pairs)
            special_base_pair_set = special_base_pairs{q};
            hex_color = convert_rgb_to_hexadecimal(special_colors{q});
            for k = 1:size(special_base_pair_set, 1);
                fprintf(fid, '(%d,%d):thickness=3,color=#%6s;', special_base_pair_set(k, 1), special_base_pair_set(k, 2), hex_color);
            end;
        end;
        fprintf(fid, '" />\n');
        
    end;
    
    if (exist('offset', 'var') || exist('bpp_values', 'var'));
        
        if exist('offset', 'var');
            fprintf(fid, '<param name="baseNum" value="#FFFFFF" />\n');
            fprintf(fid, '<param name="periodNum" value="1000" />\n');
        end;
        
        fprintf(fid, '<param name="annotations" value="');
        if exist('offset', 'var');
            PERIOD = 50;
            for i = 1:length(sequence)
                if (mod(i+offset, PERIOD) == 0)
                    fprintf(fid, '%d:type=B,anchor=%d,color=#000000,size=8;', i + offset, i);
                end;
            end;
        end;
        
        if exist('bpp_values', 'var') && ~isempty(bpp_values);
            if length(bpp_values) ~= length(bpp_anchor_bases);
                fprintf('Must specify a bpp_anchor_base for each bpp_value \n');
            end;
            for i = 1:length(bpp_values)
                %fprintf( fid, '%3.0f%%:type=L,anchor=%d,color=#303030,size=9;', 100*bpp_values(i), bpp_anchor_bases(i) );
                %fprintf( fid, '%3.0f%%:type=L,anchor=%d,color=#FF3030,size=9;', 100*bpp_values(i), bpp_anchor_bases(i) );
                fprintf(fid, '%3.0f%%:type=L,anchor=%d,color=#009000,size=9;', 100 * bpp_values(i), bpp_anchor_bases(i));
            end;
        end;
        
        fprintf(fid, '">\n');
    end;
    
    
    fprintf(fid, '%s\n', '<PARAM name="flat" value="true" />');
    fprintf(fid, '%s\n', '</applet></BODY></HTML>');

    fclose(fid);
    
    fprintf( '\n\nYou probably will not be able to run the HTML version of VARNA -- no Web browser supports it anymore!\nInstead, rerun with .html in the name of your file, and you will get an interactive version of VARNA.\nOr use .jpg or .png as the file extension.\n' );
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hex_string = convert_rgb_to_hexadecimal(rgb)

hex_string = '';
for i = 1:3;
    hex_string = [hex_string, pad_with_zero(dec2hex(floor(rgb(i)*255)))];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = pad_with_zero(s)
if length(s) == 1;
    s = ['0', s];
end;

