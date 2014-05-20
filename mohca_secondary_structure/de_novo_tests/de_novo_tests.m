% De novo tests

rdatdir = 'RDATS/';


% No pseudoknots 
outdir = 'results_no_pk/'

% HoxA9D
rdatname = '20140417_HoxA9Dnew_REBALANCE.COHCOA.rdat';
name = 'HoxA9D_rebalanced';
mohca_RNAstructure(name, [rdatdir, rdatname], [], outdir, 0, 0, 1);

rdatname = '20140417_HoxA9DnewpAadaptCprime.RAW.2.COHCOA.rdat';
name = 'HoxA9D';
mohca_RNAstructure(name, [rdatdir, rdatname], [], outdir, 0, 0, 0.5);

% truncated domain

rdatname = '20140214_Hox189.RAW.2.COHCOA.SQR.rdat';
name = 'HoxA9D_truncated';
mohca_RNAstructure(name, [rdatdir, rdatname], [], outdir, 0, 0, 1);

% With pseudoknots (SHAPEKnots)

outdir = 'results/'
% HoxA9D
rdatname = '20140417_HoxA9Dnew_REBALANCE.COHCOA.rdat';
name = 'HoxA9D_rebalanced';
mohca_RNAstructure(name, [rdatdir, rdatname], [], outdir, 1, 0, 1);

rdatname = '20140417_HoxA9DnewpAadaptCprime.RAW.2.COHCOA.rdat';
name = 'HoxA9D';
mohca_RNAstructure(name, [rdatdir, rdatname], [], outdir, 1, 0, 0.5);

% truncated domain

rdatname = '20140214_Hox189.RAW.2.COHCOA.SQR.rdat';
name = 'HoxA9D_truncated';
mohca_RNAstructure(name, [rdatdir, rdatname], [], outdir, 1, 0, 1);


