Test training
  $ python -m tcr_pmhc_predicted_structures.apps.train_pae_aggregator \
  > --log-level info \
  > --seed 123 \
  > --num-epochs 1 \
  > --batch-size 1 \
  > --eval-batch-size 1 \
  > -o models \
  > $TESTDIR/data/data.h5
  * - INFO: Setting random seed to 123 (glob)
  * - INFO: Parameter: training_data='*/data/data.h5' (glob)
  * - INFO: Parameter: output='models' (glob)
  * - INFO: Parameter: seed=123 (glob)
  * - INFO: Parameter: label_name='label' (glob)
  * - INFO: Parameter: num_epochs=1 (glob)
  * - INFO: Parameter: batch_size=1 (glob)
  * - INFO: Parameter: learning_rate=0.001 (glob)
  * - INFO: Parameter: eval_batch_size=1 (glob)
  * - INFO: Parameter: eval_interval=50 (glob)
  * - INFO: Parameter: log_level='info' (glob)
  * - INFO: Parameter: log_file=None (glob)
  * - INFO: Loading training data from */data/data.h5 (glob)
  * - INFO: Creating output directory: models (glob)
  * - INFO: Starting fold 1 (glob)
  * - INFO: Initialising model (glob)
  * - INFO: Starting training (glob)
  * - INFO: Starting epoch 1 (glob)
  * - INFO: Evaluating at iteration 1 (glob)
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  * - INFO: Loss: * (glob)
  * - INFO: Training ROC-AUC: nan (glob)
  * - INFO: Validation ROC-AUC: nan (glob)
  * - INFO: Saving model to models/model_1.pt (glob)
  * - INFO: Starting fold 2 (glob)
  * - INFO: Initialising model (glob)
  * - INFO: Starting training (glob)
  * - INFO: Starting epoch 1 (glob)
  * - INFO: Evaluating at iteration 1 (glob)
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn( (glob)
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn( (glob)
  * - INFO: Loss: * (glob)
  * - INFO: Training ROC-AUC: nan (glob)
  * - INFO: Validation ROC-AUC: nan (glob)
  * - INFO: Saving model to models/model_2.pt (glob)
  * - INFO: Starting fold 3 (glob)
  * - INFO: Initialising model (glob)
  * - INFO: Starting training (glob)
  * - INFO: Starting epoch 1 (glob)
  * - INFO: Evaluating at iteration 1 (glob)
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  *: UndefinedMetricWarning: Only one class is present in y_true. ROC AUC score is not defined in that case. (glob)
    warnings.warn(
  * - INFO: Loss: * (glob)
  * - INFO: Training ROC-AUC: nan (glob)
  * - INFO: Validation ROC-AUC: nan (glob)
  * - INFO: Saving model to models/model_3.pt (glob)

  $ python <<EOF
  > import os
  > import h5py
  > import numpy as np
  > import torch
  > from tcr_pmhc_predicted_structures.models import PAEAggregator
  > with h5py.File(os.path.join(os.environ['TESTDIR'], 'data', 'data.h5')) as fh:
  >     paes = fh['pae'][:]
  > 
  > ref_model1 = PAEAggregator(3120)
  > ref_model1.load_state_dict(torch.load(
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'models', 'model_1.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > ref_model1.eval()
  > ref_predictions1 = ref_model1(torch.tensor(paes, dtype=torch.float32))
  > test_model1 = PAEAggregator(3120)
  > test_model1.load_state_dict(torch.load(
  >     os.path.join('models', 'model_1.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > test_model1.eval()
  > test_predictions1 = test_model1(torch.tensor(paes, dtype=torch.float32))
  > np.testing.assert_array_almost_equal(test_predictions1.detach().numpy(), ref_predictions1.detach().numpy())
  > ref_model2 = PAEAggregator(3120)
  > ref_model2.load_state_dict(torch.load(
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'models', 'model_2.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > ref_model2.eval()
  > ref_predictions2 = ref_model2(torch.tensor(paes, dtype=torch.float32))
  > test_model2 = PAEAggregator(3120)
  > test_model2.load_state_dict(torch.load(
  >     os.path.join('models', 'model_2.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > test_model2.eval()
  > test_predictions2 = test_model2(torch.tensor(paes, dtype=torch.float32))
  > np.testing.assert_array_almost_equal(test_predictions2.detach().numpy(), ref_predictions2.detach().numpy())
  > ref_model3 = PAEAggregator(3120)
  > ref_model3.load_state_dict(torch.load(
  >     os.path.join(os.environ['TESTDIR'], 'reference', 'models', 'model_3.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > ref_model3.eval()
  > ref_predictions3 = ref_model3(torch.tensor(paes, dtype=torch.float32))
  > test_model3 = PAEAggregator(3120)
  > test_model3.load_state_dict(torch.load(
  >     os.path.join('models', 'model_3.pt'),
  >     map_location=torch.device('cpu'),
  >     weights_only=True,
  > ))
  > test_model3.eval()
  > test_predictions3 = test_model3(torch.tensor(paes, dtype=torch.float32))
  > np.testing.assert_array_almost_equal(test_predictions3.detach().numpy(), ref_predictions3.detach().numpy())
  > EOF
