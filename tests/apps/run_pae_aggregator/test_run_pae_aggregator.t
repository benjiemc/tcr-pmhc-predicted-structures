Test untrained.
  $ python -m tcr_pmhc_predicted_structures.apps.run_pae_aggregator \
  > --mode untrained \
  > --paes $TESTDIR/data/paes.h5 \
  > --pcas $TESTDIR/data/pcas.csv \
  > --model $TESTDIR/data/model \
  > -o untrained.csv

  $ python <<EOF
  > import os
  > import pandas as pd
  > predictions = pd.read_csv('untrained.csv')
  > reference = pd.read_csv(os.path.join(os.environ['TESTDIR'], 'reference', 'untrained.csv'))
  > pd.testing.assert_frame_equal(predictions, reference, check_exact=False)
  > EOF

Test evaluation.
  $ python -m tcr_pmhc_predicted_structures.apps.run_pae_aggregator \
  > --mode evaluation \
  > --paes $TESTDIR/data/paes.h5 \
  > --pcas $TESTDIR/data/pcas.csv \
  > --model $TESTDIR/data/model \
  > -o evaluation.csv

  $ python <<EOF
  > import os
  > import pandas as pd
  > predictions = pd.read_csv('evaluation.csv')
  > reference = pd.read_csv(os.path.join(os.environ['TESTDIR'], 'reference', 'evaluation.csv'))
  > pd.testing.assert_frame_equal(predictions, reference, check_exact=False)
  > EOF

Test run.
  $ python -m tcr_pmhc_predicted_structures.apps.run_pae_aggregator \
  > --mode run \
  > --paes $TESTDIR/data/paes.h5 \
  > --pcas $TESTDIR/data/pcas.csv \
  > --model $TESTDIR/data/model \
  > -o run.csv

  $ python <<EOF
  > import os
  > import pandas as pd
  > predictions = pd.read_csv('run.csv')
  > reference = pd.read_csv(os.path.join(os.environ['TESTDIR'], 'reference', 'run.csv'))
  > pd.testing.assert_frame_equal(predictions, reference, check_exact=False)
  > EOF
