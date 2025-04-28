Test.
  $ python -m tcr_pmhc_predicted_structures.apps.create_af3_json \
  > --log-level error \
  > --mhc-sequences $TESTDIR/data/hla_a.json $TESTDIR/data/h2_d.json \
  > -o jsons \
  > $TESTDIR/data/sequences.csv

  $ jq --sort-keys . jsons/0.json > test
  $ jq --sort-keys . $TESTDIR/reference/0.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/1.json > test
  $ jq --sort-keys . $TESTDIR/reference/1.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/2.json > test
  $ jq --sort-keys . $TESTDIR/reference/2.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/3.json > test
  $ jq --sort-keys . $TESTDIR/reference/3.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/4.json > test
  $ jq --sort-keys . $TESTDIR/reference/4.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/5.json > test
  $ jq --sort-keys . $TESTDIR/reference/5.json > ref
  $ diff test ref
