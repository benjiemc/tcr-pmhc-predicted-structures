Test app
  $ python -m tcr_pmhc_predicted_structures.apps.create_contact_map \
  > -o contacts.csv \
  > $TESTDIR/data/tcr_pmhc.pdb

  $ diff contacts.csv $TESTDIR/reference/contacts.csv
