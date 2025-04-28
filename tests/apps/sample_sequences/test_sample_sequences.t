Test sampling.
  $ python -m tcr_pmhc_predicted_structures.apps.sample_sequences \
  > --seed 123 \
  > --sample-size 5 \
  > -o sample.csv \
  > $TESTDIR/data/sequences.csv

  $ diff sample.csv $TESTDIR/reference/sample.csv

Test sampling and background generation.
  $ python -m tcr_pmhc_predicted_structures.apps.sample_sequences \
  > --seed 123 \
  > --generate-background \
  > --background-proportion 1 \
  > --sample-size 5 \
  > -o sample_with_background.csv \
  > $TESTDIR/data/sequences.csv

  $ diff sample_with_background.csv $TESTDIR/reference/sample_with_background.csv
