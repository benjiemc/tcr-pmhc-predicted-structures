Test app (standard)
  $ python -m tcr_pmhc_predicted_structures.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --output-summary-csv structure_summary.csv \
  > -o output \
  > $TESTDIR/data/stcrdab-mock/

  $ diff structure_summary.csv $TESTDIR/reference/standard/structure_summary.csv

  $ diff output/3qiw_CDEAB.pdb $TESTDIR/reference/standard/3qiw_CDEAB.pdb
  $ diff output/7q9b_DECAB.pdb $TESTDIR/reference/standard/7q9b_DECAB.pdb
  $ diff output/7q9b_IJHFG.pdb $TESTDIR/reference/standard/7q9b_IJHFG.pdb

Test removing structures missing residues
  $ python -m tcr_pmhc_predicted_structures.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --remove-structures-missing-residues \
  > --output-summary-csv no_missing_resiudes_structure_summary.csv \
  > -o output-no-missing-residues \
  > $TESTDIR/data/stcrdab-mock/

  $ diff no_missing_resiudes_structure_summary.csv $TESTDIR/reference/output-no-missing-residues/structure_summary.csv

Test fixing structures missing reisdues
  $ python -m tcr_pmhc_predicted_structures.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --remove-structures-missing-residues \
  > --fix-structures-missing-residues \
  > --output-summary-csv fix_missing_residues_structure_summary.csv \
  > -o output-fixed-missing-residues \
  > $TESTDIR/data/stcrdab-mock/

  $ diff fix_missing_residues_structure_summary.csv $TESTDIR/reference/fixed-missing-residues/structure_summary.csv
  $ diff output-fixed-missing-residues/3qiw_CDEAB.pdb $TESTDIR/reference/fixed-missing-residues/3qiw_CDEAB.pdb
  $ diff output-fixed-missing-residues/6v19_DECAB.pdb $TESTDIR/reference/fixed-missing-residues/6v19_DECAB.pdb

Test structural similarity cutoff
  $ python -m tcr_pmhc_predicted_structures.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --structural-similarity-cutoff 2.0 \
  > --output-summary-csv structural_similarity_cutoff_structure_summary.csv \
  > -o output-structural-similarity-cutoff \
  > $TESTDIR/data/stcrdab-mock/

  $ diff structural_similarity_cutoff_structure_summary.csv $TESTDIR/reference/structural-similarity-cutoff/structure_summary.csv

Test crop structures
  $ python -m tcr_pmhc_predicted_structures.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --crop-structures \
  > --output-summary-csv crop_structure_structure_summary.csv \
  > -o output-crop \
  > $TESTDIR/data/stcrdab-mock/

  $ diff crop_structure_structure_summary.csv $TESTDIR/reference/crop/structure_summary.csv

  $ diff output-crop/3qiw_CDEAB.pdb $TESTDIR/reference/crop/3qiw_CDEAB.pdb
  $ diff output-crop/7q9b_DECAB.pdb $TESTDIR/reference/crop/7q9b_DECAB.pdb
  $ diff output-crop/7q9b_IJHFG.pdb $TESTDIR/reference/crop/7q9b_IJHFG.pdb

Test remove HETATMs
  $ python -m tcr_pmhc_predicted_structures.apps.select_stcrdab_tcr_pmhc_structures \
  > --log-level error \
  > --tcr-types abTCR \
  > --mhc-types MH1 MH2 \
  > --antigen-types peptide \
  > --remove-het-atoms \
  > --output-summary-csv remove_hetatms_structure_summary.csv \
  > -o output-no-het-atoms \
  > $TESTDIR/data/stcrdab-mock/

  $ diff remove_hetatms_structure_summary.csv $TESTDIR/reference/no-het-atoms/structure_summary.csv

  $ diff output-no-het-atoms/3qiw_CDEAB.pdb $TESTDIR/reference/no-het-atoms/3qiw_CDEAB.pdb
  $ diff output-no-het-atoms/7q9b_DECAB.pdb $TESTDIR/reference/no-het-atoms/7q9b_DECAB.pdb
  $ diff output-no-het-atoms/7q9b_IJHFG.pdb $TESTDIR/reference/no-het-atoms/7q9b_IJHFG.pdb
