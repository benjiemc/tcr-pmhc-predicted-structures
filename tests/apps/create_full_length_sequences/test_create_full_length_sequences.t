Test app.
  $ python -m tcr_pmhc_predicted_structures.apps.create_full_length_sequences \
  > --mhc-sequences $TESTDIR/data/hla_a.json $TESTDIR/data/h2_d.json \
  > -o test.csv \
  > $TESTDIR/data/sequences.csv

  $ diff test.csv $TESTDIR/reference/sequences.csv

Test skipping unstichable row.
  $ python -m tcr_pmhc_predicted_structures.apps.create_full_length_sequences \
  > --mhc-sequences $TESTDIR/data/hla_a.json \
  > -o skipped.csv \
  > $TESTDIR/data/sequences_unstitchable.csv
  * - WARNING: Unable to stitch: *. Skipping... (glob)

  $ cut -d, -f1 skipped.csv
  id
  good
