  $ python -m tcr_pmhc_predicted_structures.apps.create_pae_dataset \
  > --sequences $TESTDIR/data/sequences.csv \
  > --confidences $TESTDIR/data/10860_1_0_pae.txt  $TESTDIR/data/10860_1_1_pae.txt  $TESTDIR/data/4mxq_CD_7byd_CAB_1_0_pae.txt  $TESTDIR/data/66254_1_4_pae.txt \
  > -o dataset.h5

  $ python <<EOF
  > import os
  > 
  > import numpy as np
  > import h5py
  > 
  > with h5py.File('dataset.h5', 'r') as fh:
  >     ids = np.array([name.decode() for name in fh['id'][:]])
  >     seeds = fh['seed'][:]
  >     samples = fh['sample'][:]
  >     labels = fh['label'][:]
  >     sources = fh['source'][:]
  >     folds = fh['fold'][:]
  >     paes = fh['pae'][:]
  > 
  > ref_ids = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'id.txt'), dtype=str)
  > ref_seeds = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'seed.txt'), dtype=int)
  > ref_samples = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'sample.txt'), dtype=int)
  > ref_labels = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'label.txt'), dtype=int)
  > ref_sources = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'source.txt'), dtype=int)
  > ref_folds = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'fold.txt'), dtype=int)
  > ref_paes = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pae.txt'))
  > 
  > np.testing.assert_array_equal(ids, ref_ids)
  > np.testing.assert_array_equal(seeds, ref_seeds)
  > np.testing.assert_array_equal(samples, ref_samples)
  > np.testing.assert_array_equal(labels, ref_labels)
  > np.testing.assert_array_equal(sources, ref_sources)
  > np.testing.assert_array_equal(folds, ref_folds)
  > np.testing.assert_array_almost_equal(paes, ref_paes)
  > EOF

Test json format
  $ python -m tcr_pmhc_predicted_structures.apps.create_pae_dataset \
  > --sequences $TESTDIR/data/sequences.csv \
  > --format json \
  > --confidences $TESTDIR/data/10860_seed-1_sample-0_confidences.json  $TESTDIR/data/10860_seed-1_sample-1_confidences.json  $TESTDIR/data/4mxq_CD_7byd_CAB_seed-1_sample-0_confidences.json  $TESTDIR/data/66254_seed-1_sample-4_confidences.json \
  > -o dataset_from_json.h5

  $ python <<EOF
  > import os
  > 
  > import numpy as np
  > import h5py
  > 
  > with h5py.File('dataset_from_json.h5', 'r') as fh:
  >     ids = np.array([name.decode() for name in fh['id'][:]])
  >     seeds = fh['seed'][:]
  >     samples = fh['sample'][:]
  >     labels = fh['label'][:]
  >     sources = fh['source'][:]
  >     folds = fh['fold'][:]
  >     paes = fh['pae'][:]
  > 
  > ref_ids = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'id.txt'), dtype=str)
  > ref_seeds = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'seed.txt'), dtype=int)
  > ref_samples = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'sample.txt'), dtype=int)
  > ref_labels = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'label.txt'), dtype=int)
  > ref_sources = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'source.txt'), dtype=int)
  > ref_folds = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'fold.txt'), dtype=int)
  > ref_paes = np.loadtxt(os.path.join(os.environ['TESTDIR'], 'reference', 'pae.txt'))
  > 
  > np.testing.assert_array_equal(ids, ref_ids)
  > np.testing.assert_array_equal(seeds, ref_seeds)
  > np.testing.assert_array_equal(samples, ref_samples)
  > np.testing.assert_array_equal(labels, ref_labels)
  > np.testing.assert_array_equal(sources, ref_sources)
  > np.testing.assert_array_equal(folds, ref_folds)
  > np.testing.assert_array_almost_equal(paes, ref_paes)
  > EOF