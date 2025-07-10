  $ python -m tcr_pmhc_predicted_structures.apps.extract_confidence_from_af3 \
  > -o test.csv \
  > --summary-confidence-paths $TESTDIR/data/4586_seed-1_sample-0_summary_confidences.json $TESTDIR/data/4586_seed-1_sample-1_summary_confidences.json $TESTDIR/data/4586_seed-1_sample-2_summary_confidences.json $TESTDIR/data/4586_seed-1_sample-3_summary_confidences.json $TESTDIR/data/4586_seed-1_sample-4_summary_confidences.json \
  > --confidence-paths $TESTDIR/data/4586_seed-1_sample-0_confidences.json $TESTDIR/data/4586_seed-1_sample-1_confidences.json $TESTDIR/data/4586_seed-1_sample-2_confidences.json $TESTDIR/data/4586_seed-1_sample-3_confidences.json $TESTDIR/data/4586_seed-1_sample-4_confidences.json \
  > --sequence-paths $TESTDIR/data/4586.json
  Error: Something really went wrong that has not been handled
  'C'
  Error: Something really went wrong that has not been handled
  'C'
  Error: Something really went wrong that has not been handled
  'C'
  Error: Something really went wrong that has not been handled
  'C'
  Error: Something really went wrong that has not been handled
  'C'
  Error: Something really went wrong that has not been handled
  'C'
  Error: Something really went wrong that has not been handled
  'C'
  Error: Something really went wrong that has not been handled
  'C'
  Error: Something really went wrong that has not been handled
  'C'
  Error: Something really went wrong that has not been handled
  'C'

  $ diff test.csv $TESTDIR/reference/confidence_scores.csv
