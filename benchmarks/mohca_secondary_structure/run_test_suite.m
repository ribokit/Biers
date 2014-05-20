function run_test_suite(name, rdatname, structure, outdir)
threshold = 1.2;

% No data
run_test_case(name, rdatname, structure, [outdir, 'no_data/'], -1, false, false, false);
% Data as is
% Just diagonal filter
run_test_case(name, rdatname, structure, [outdir, 'diagonal/'], threshold, false, false, false);

% Diagonal filter + combine diagonals
run_test_case(name, rdatname, structure, [outdir, 'diagonal_combine/'], threshold, true, false, false);

% Diagonal filter + combine diagonals + broaden
run_test_case(name, rdatname, structure, [outdir, 'diagonal_combine_broaden/'], threshold, true, true, false);

% Diagonal filter + broaden
run_test_case(name, rdatname, structure, [outdir, 'diagonal_combine_broaden/'], threshold, false, true, false);

% Binarized data
% Just diagonal filter
run_test_case(name, rdatname, structure, [outdir, 'diagonal_binarized/'], threshold, false, false, true);

% Diagonal filter + combine diagonals
run_test_case(name, rdatname, structure, [outdir, 'diagonal_combine_binarized/'], threshold, true, false, true);

% Diagonal filter + combine diagonals + broaden
run_test_case(name, rdatname, structure, [outdir, 'diagonal_combine_broaden_binarized/'], threshold, true, true, true);

% Diagonal filter + broaden
run_test_case(name, rdatname, structure, [outdir, 'diagonal_combine_broaden/'], threshold, false, true, true);
end