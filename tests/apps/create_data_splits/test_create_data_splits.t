Test splitting on pMHCs
  $ python -m tcr_pmhc_predicted_structures.apps.create_data_splits \
  > --seed 123 \
  > --split-type pMHC \
  > -o sequences_pmhc_split.csv \
  > $TESTDIR/data/sequences.csv
  * - WARNING: Insufficient data to create 5 folds. Only 3 fold(s) will be created (glob)

  $ diff sequences_pmhc_split.csv $TESTDIR/reference/sequences_pmhc_split.csv

Test splitting on TCRs
  $ python -m tcr_pmhc_predicted_structures.apps.create_data_splits \
  > --seed 123 \
  > --split-type tcr \
  > -o sequences_tcr_split.csv \
  > $TESTDIR/data/sequences.csv

  $ diff sequences_tcr_split.csv $TESTDIR/reference/sequences_tcr_split.csv

Test split on levenshtein distance
  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > with h5py.File('distances.h5', 'w') as hdf5:
  >     for entity in [
  >         'cdr1_alpha', 'cdr2_alpha', 'cdr3_alpha', 'cdr1_beta', 'cdr2_beta', 'cdr3_beta', 'peptide', 'mhc_pseudo'
  >     ]:
  >         group = hdf5.create_group(entity)
  >         with open(os.path.join(os.environ['TESTDIR'], 'data', f'{entity}s.txt')) as fh:
  >             group['names'] = [line.strip() for line in fh.readlines() if line]
  >         group['distance_matrix'] = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'data', f'{entity}_distances.txt'), dtype=int)
  > EOF

  $ python -m tcr_pmhc_predicted_structures.apps.create_data_splits \
  > --seed 123 \
  > --split-type levenshtein \
  > --split-entities cdr1_alpha cdr3_alpha \
  > --split-distance 5 \
  > --distances distances.h5 \
  > -o sequences_levenshtein_split.csv \
  > $TESTDIR/data/sequences.csv

  $ diff sequences_levenshtein_split.csv $TESTDIR/reference/sequences_levenshtein_split.csv
