Test app.
  $ python -m tcr_pmhc_predicted_structures.apps.collate_sequence_data \
  > --log-level error \
  > --iedb-path $TESTDIR/data/iedb.csv \
  > --vdjdb-path $TESTDIR/data/vdjdb.tsv \
  > --itrap-path $TESTDIR/data/itrap.csv \
  > --mcpas-tcr-path $TESTDIR/data/mcpas_tcr.csv \
  > --mhc-sequences $TESTDIR/data/h2_d.json $TESTDIR/data/hla_a.json $TESTDIR/data/hla_b.json \
  > -o sequences.csv

  $ diff sequences.csv $TESTDIR/reference/sequences.csv
