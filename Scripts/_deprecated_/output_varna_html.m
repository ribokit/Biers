function [special_base_pairs, special_colors] = output_varna_html(varna_file, sequence, display_structure, native_structure, modeled_structure, offset, mutpos, crystpos, shape, bpp, color_mode)
% [special_base_pairs, special_colors] = output_varna( ...
%                               varna_file, sequence, display_structure, native_structure, modeled_structure, 
%                               offset, mutpos, crystpos, shape, bpp, color_mode); 
%
% Generate HTML file with VARNA applet for RNA secstr visualziation.
% **Make sure to set your path in get_varna_exe.m
%
% [Input]
% varna_file            Required            File name (must end in: .html)                                                          .png
% sequence              Required            RNA sequence for display
% display_structure     Required            RNA secstr for display
% native_structure      Optional            Reference/native secstr used for secstr comparison.
% modeled_structure     Optional            Predicted/new secstr used for secstr comparison.
% offset                Optional            Sequence numbering offset, default 0.
% mutpos                Optional            [?]
% crystpos              Optional            [?]
% shape                 Optional            SHAPE reactivity profile for nucleotide coloring, default none.
% bpp                   Optional            Base-pairing probability matrix for helix-wise confidence score.
% color_mode            Optional            Secstr comparison line color set, default 1.
%
% [Output]
% special_base_pairs    Index of base-pairs that differ between native_structure and modeled_structure
% special_colors        Secstr comparison color set
%
% by T47, 2013-2015
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist('native_structure','var') || isempty(native_structure); native_structure = display_structure; end;
if ~exist('modeled_structure','var') || isempty(modeled_structure); modeled_structure = display_structure; end;
if ~exist('bpp', 'var') || isempty(bpp); bpp = []; end;
if ~exist('shape', 'var') || isempty(shape); shape = []; end;

if ~exist('offset','var') || isempty(offset); offset = 0; end;
if ~exist('mutpos','var') || isempty(mutpos); mutpos = [1:length(sequence)] + offset; end;
if ~exist('crystpos','var') || isempty(crystpos); crystpos = [1:length(sequence)] + offset; end;
if ~exist('color_mode', 'var') || isempty(color_mode); color_mode = 1; end;

output_varna( varna_file, sequence, display_structure, native_structure, modeled_structure, offset, mutpos, crystpos, shape, bpp, color_mode );

fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n' );
fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n' );
fprintf( [mfilename, ' is deprecated. Use output_varna() instead.\n']  );
fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n' );
fprintf( 'WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n' );


