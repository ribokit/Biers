% Clarence Cheng
% May 17, 2014
% Develop MOHCA-seq data extraction script to narrow down data for conversion into pseudoenergies 

% Prepare data and 
r_p4p6      = read_rdat_file('~/Dropbox/MOHCAseqPAPER/Data/ALLRDATS/P4P6.COMBINED.COHCOA.SQR.rdat');
D_p4p6      = r_p4p6.reactivity;
D_err_p4p6  = r_p4p6.reactivity_error;
ligpos_p4p6 = get_ligpos(r_p4p6);
seqpos_p4p6 = r_p4p6.seqpos;
secstr_p4p6 = {'GGCCAAAACAACGGAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCAAAACCAAACCAAAGAAACAACAACAACAACX', ...
               '.................((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...))..................................'};
mohcaplot(D_p4p6, seqpos_p4p6, ligpos_p4p6, 'test', '', '', secstr_p4p6);

r_add = read_rdat_file('~/Dropbox/MOHCAseqPAPER/Data/5_ADD_RBSW/A1_Lig_Asc/COMBINED.COHCOA.SQR.rdat');
D_add = r_add.reactivity;
D_err_add = r_add.reactivity_error;
ligpos_add = get_ligpos(r_add);
seqpos_add = r_add.seqpos;
secstr_add = {'GGAAAGCAAUUCGAGUAGAAUUGGAAAGGGAAAGAAACGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUGAAAACAAAGUUAAGGAGUACUUAACACAAAGAAACAACAACAACAACX', ...
              '......((((((.....))))))..............(((((((((...((((((.........))))))........((((((.......))))))..)))))))))........((((((.....)))))).......................'};
mohcaplot(D_add, seqpos_add, ligpos_add, 'test', '','', secstr_add);

r_fngly      = read_rdat_file('~/Dropbox/MOHCAseqPAPER/Data/7_FN_RBSW/A1_Lig_Asc/COMBINED.COHCOA.SQR.rdat');
D_fngly      = r_fngly.reactivity;
D_err_fngly  = r_fngly.reactivity_error;
ligpos_fng   = get_ligpos(r_fngly);
seqpos_fng   = r_fngly.seqpos;
secstr_fng = { 'GGCAAUUCGAGUAGAAUUGACAGAGAGGAUAUGAGGAGAGAUUUCAUUUUAAUGAAACACCGAAGAAGUAAAUCUUUCAGGUAAAAAGGACUCAUAUUGGACGAACCUCUGGAGAGCUUAUCUAAGAGAUAACACCGAAGGAGCAAAGCUAAUUUUAGCCUAAACUCUCAGGUAAAAGGACGGAGAAAACACAAGUUCAGGAGUACUGAACCAAAGAAACAACAACAACAACX', ...
               '..((((((.....))))))........((((((((......((((((....)))))).(((...((((.....))))..)))........))))))))........(((((......(((((.....))))).(((...((((....((((....)))).....))))..))).......))))).........((((((.....))))))......................'};
mohcaplot(D_fngly, seqpos_fng, ligpos_fng, 'test', '','', secstr_fng);

%% P4-P6
% Vertical comparison
D_comb_p4p6 = combine_diag(D_p4p6, 1, [0.1 0.25 0.5 0.75 1], 7, seqpos_p4p6, ligpos_p4p6, secstr_p4p6);
D_comb_add = combine_diag(D_add, 1, [0.1 0.25 0.5 0.75 1], 7, seqpos_add, ligpos_add, secstr_add);
D_comb_fng = combine_diag(D_fngly, 1, [0.1 0.25 0.5 0.75 1], 7, seqpos_fng, ligpos_fng, secstr_fng);


D_test = mohca_symmetric(D_p4p6);
D_filt = diagonal_filter(D_test);
D_comb = combine_diag(D_filt, 1, [0.1 0.25 0.5 0.75 1], 7, seqpos_p4p6, ligpos_p4p6, secstr_p4p6);

D_comb = combine_diag(D_filt, 1, 0.025, 7, seqpos_p4p6, ligpos_p4p6, secstr_p4p6);

% Try diagonal comparison
D_comb = combine_diag(D_filt, 1, 0.025, 2, seqpos_p4p6, ligpos_p4p6, secstr_p4p6);
D_comb2 = combine_diag(D_filt, '', 0.025, 2, seqpos_p4p6, ligpos_p4p6, secstr_p4p6);
D_comb_comb = combine_diag(D_comb, '', 0.025, 2, seqpos_p4p6, ligpos_p4p6, secstr_p4p6);

D_comb = combine_diag(D_filt, 1, 0.025, 2, seqpos_p4p6, ligpos_p4p6, secstr_p4p6);
D_comb_comb = combine_diag(D_comb, '', 0.025, 2, seqpos_p4p6, ligpos_p4p6, secstr_p4p6);

%% FN glycine riboswitch
% Final comparison for the day
D_sym_fng = mohca_symmetric(D_fngly);
D_filt_fng = diagonal_filter(D_sym_fng);
D_comb_fng = combine_diag(D_filt_fng, 1, 0.025, 2, seqpos_fng, ligpos_fng, secstr_fng);
D_comb_comb_fng = combine_diag(D_comb_fng, '', 0.025, 2, seqpos_fng, ligpos_fng, secstr_fng);


