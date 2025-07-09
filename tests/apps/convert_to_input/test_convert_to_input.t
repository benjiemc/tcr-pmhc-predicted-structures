Create AF3 output
  $ python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
  > --log-level error \
  > --format alphafold3 \
  > -o jsons \
  > $TESTDIR/data/sequences.csv

  $ jq --sort-keys . jsons/0.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/0.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/1.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/1.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/2.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/2.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/3.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/3.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/4.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/4.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons/5.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/5.json > ref
  $ diff test ref

Create AF3 inputs with paths
  $ python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
  > --log-level error \
  > --format alphafold3 \
  > --sequence-key $TESTDIR/data/sequence_key.txt \
  > --relative-paths \
  > --unpaired-msa-paths $TESTDIR/data/msas/*_unpaired.a3m \
  > --paired-msa-paths $TESTDIR/data/msas/*_paired.a3m \
  > --template-paths $TESTDIR/data/templates/* \
  > -o jsons-paths \
  > $TESTDIR/data/sequences.csv

  $ jq --sort-keys . jsons-paths/0.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/0_path.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons-paths/1.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/1_path.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons-paths/2.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/2_path.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons-paths/3.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/3_path.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons-paths/4.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/4_path.json > ref
  $ diff test ref

  $ jq --sort-keys . jsons-paths/5.json > test
  $ jq --sort-keys . $TESTDIR/reference/alphafold3/5_path.json > ref
  $ diff test ref

Create chai input
  $ python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
  > --format chai \
  > -o chai \
  > $TESTDIR/data/sequences.csv

  $ diff chai/0.fasta $TESTDIR/reference/chai/0.fasta
  $ diff chai/1.fasta $TESTDIR/reference/chai/1.fasta
  $ diff chai/2.fasta $TESTDIR/reference/chai/2.fasta
  $ diff chai/3.fasta $TESTDIR/reference/chai/3.fasta
  $ diff chai/4.fasta $TESTDIR/reference/chai/4.fasta
  $ diff chai/5.fasta $TESTDIR/reference/chai/5.fasta

Create chai input with restraints.
  $ python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
  > --format chai \
  > --chai-restraints $TESTDIR/data/restraints.csv \
  > -o chai-with-restraints \
  > $TESTDIR/data/sequences.csv

  $ diff chai-with-restraints/0.fasta $TESTDIR/reference/chai/0.fasta
  $ diff chai-with-restraints/1.fasta $TESTDIR/reference/chai/1.fasta
  $ diff chai-with-restraints/2.fasta $TESTDIR/reference/chai/2.fasta
  $ diff chai-with-restraints/3.fasta $TESTDIR/reference/chai/3.fasta
  $ diff chai-with-restraints/4.fasta $TESTDIR/reference/chai/4.fasta
  $ diff chai-with-restraints/5.fasta $TESTDIR/reference/chai/5.fasta

  $ diff chai-with-restraints/0_restraints.csv $TESTDIR/reference/chai/0_restraints.csv
  $ diff chai-with-restraints/1_restraints.csv $TESTDIR/reference/chai/1_restraints.csv
  $ diff chai-with-restraints/2_restraints.csv $TESTDIR/reference/chai/2_restraints.csv
  $ diff chai-with-restraints/3_restraints.csv $TESTDIR/reference/chai/3_restraints.csv
  $ diff chai-with-restraints/4_restraints.csv $TESTDIR/reference/chai/4_restraints.csv
  $ diff chai-with-restraints/5_restraints.csv $TESTDIR/reference/chai/5_restraints.csv
