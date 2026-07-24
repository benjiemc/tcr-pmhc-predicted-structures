Rows that stitchr cannot build are skipped with a warning rather than
aborting the whole run (regression: create_full_length_sequences raised
KeyError on the non-existent 'tcr_species' column instead of logging a skip).

  $ python -m tcr_pmhc_predicted_structures.apps.create_full_length_sequences \
  > --mhc-sequences $TESTDIR/data/hla_a.json \
  > -o out.csv \
  > $TESTDIR/data/sequences_unstitchable.csv 2> stderr.log

The unstitchable ('bad') row is skipped; the stitchable ('good') row is kept:

  $ cut -d, -f1 out.csv
  id
  good

And the skip is reported, not fatal:

  $ grep -c "Unable to stitch" stderr.log
  1
