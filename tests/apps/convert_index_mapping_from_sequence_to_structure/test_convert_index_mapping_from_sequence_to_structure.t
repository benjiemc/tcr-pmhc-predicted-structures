Test app.
  $ python -m tcr_pmhc_predicted_structures.apps.convert_index_mapping_from_sequence_to_structure \
  > $TESTDIR/data/6d7g.cif \
  > $TESTDIR/data/indices.json \
  > -o indices.json

  $ jq --sort-keys . indices.json > test
  $ jq --sort-keys . $TESTDIR/reference/indices.json > ref
  $ diff test ref
