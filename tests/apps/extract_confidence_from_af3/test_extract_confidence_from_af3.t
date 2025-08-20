  $ python -m tcr_pmhc_predicted_structures.apps.extract_confidence_from_af3 \
  > -o test.csv \
  > $TESTDIR/data/4586
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