Test extract TCR:pMHC PAE.
  $ python -m tcr_pmhc_predicted_structures.apps.extract_tcr_pmhc_confidence \
  > -o tcr_pmhc_pae.txt \
  > --score pae \
  > $TESTDIR/data/pae.txt \
  > $TESTDIR/data/sequences.json
  *: UserWarning: Region 1 of the B2M molecule was not identified by ANARCI and is therefore being numbered sequentially. See https://github.com/npqst/anarci-mhc/issues/6 for more information. (glob)
    warnings.warn(warning_message)

  $ diff tcr_pmhc_pae.txt $TESTDIR/reference/tcr_pmhc_pae.txt

Test extract TCR:pMHC PAE binding-only.
  $ python -m tcr_pmhc_predicted_structures.apps.extract_tcr_pmhc_confidence \
  > -o cdr_pmhc_pae.txt \
  > --score pae \
  > --binding-only \
  > $TESTDIR/data/pae.txt \
  > $TESTDIR/data/sequences.json
  *: UserWarning: Region 1 of the B2M molecule was not identified by ANARCI and is therefore being numbered sequentially. See https://github.com/npqst/anarci-mhc/issues/6 for more information. (glob)
    warnings.warn(warning_message)

  $ diff cdr_pmhc_pae.txt $TESTDIR/reference/cdr_pmhc_pae.txt

Test extract TCR:pMHC atom plddt.
  $ python -m tcr_pmhc_predicted_structures.apps.extract_tcr_pmhc_confidence \
  > -o tcr_pmhc_plddt.txt \
  > --score atom_plddt \
  > $TESTDIR/data/atom_plddt.txt \
  > $TESTDIR/data/sequences.json
  *: UserWarning: Region 1 of the B2M molecule was not identified by ANARCI and is therefore being numbered sequentially. See https://github.com/npqst/anarci-mhc/issues/6 for more information. (glob)
    warnings.warn(warning_message)

  $ diff tcr_pmhc_plddt.txt $TESTDIR/reference/tcr_pmhc_plddt.txt

Test extract TCR:pMHC atom plddt binding-only.
  $ python -m tcr_pmhc_predicted_structures.apps.extract_tcr_pmhc_confidence \
  > -o cdr_pmhc_plddt.txt \
  > --score atom_plddt \
  > --binding-only \
  > $TESTDIR/data/atom_plddt.txt \
  > $TESTDIR/data/sequences.json
  *: UserWarning: Region 1 of the B2M molecule was not identified by ANARCI and is therefore being numbered sequentially. See https://github.com/npqst/anarci-mhc/issues/6 for more information. (glob)
    warnings.warn(warning_message)

  $ diff cdr_pmhc_plddt.txt $TESTDIR/reference/cdr_pmhc_plddt.txt

Test padding for PAE.
  $ python -m tcr_pmhc_predicted_structures.apps.extract_tcr_pmhc_confidence \
  > -o pad_cdr_pmhc_pae.txt \
  > --score pae \
  > --pad \
  > --binding-only \
  > $TESTDIR/data/pae.txt \
  > $TESTDIR/data/sequences.json
  *: UserWarning: Region 1 of the B2M molecule was not identified by ANARCI and is therefore being numbered sequentially. See https://github.com/npqst/anarci-mhc/issues/6 for more information. (glob)
    warnings.warn(warning_message)

  $ diff pad_cdr_pmhc_pae.txt $TESTDIR/reference/pad_cdr_pmhc_pae.txt

Test padding for atom pLDDT.
  $ python -m tcr_pmhc_predicted_structures.apps.extract_tcr_pmhc_confidence \
  > -o pad_cdr_pmhc_atom_plddt.txt \
  > --score atom_plddt \
  > --pad \
  > --aggregate-atom-level-confidence \
  > --binding-only \
  > $TESTDIR/data/atom_plddt.txt \
  > $TESTDIR/data/sequences.json
  *: UserWarning: Region 1 of the B2M molecule was not identified by ANARCI and is therefore being numbered sequentially. See https://github.com/npqst/anarci-mhc/issues/6 for more information. (glob)
    warnings.warn(warning_message)

  $ diff pad_cdr_pmhc_atom_plddt.txt $TESTDIR/reference/pad_cdr_pmhc_atom_plddt.txt
