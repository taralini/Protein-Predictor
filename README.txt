
Protein Expression & Purification Predictor
==========================================

Contents
--------
- predictor.py : Python module exposing `ProteinExpressionPredictor`.
- clf_expr.pkl, tag_n.pkl, tag_c.pkl : Trained models.
- method_models.pkl (may be empty) and method_success_rates.json : Fallback rates.
- reg_*.pkl : Optional regressors for buffer pH, salt, and yield (trained if data available).
- run_predict.py : CLI to predict from a sequence or FASTA file.

Quick start
-----------
1) Download this folder.
2) In a terminal: `python3 run_predict.py MHHHHHHSSGVDLGTENLYFQSMAS...`  (or provide a FASTA path)
3) The tool prints JSON with:
   - expression_success_prob
   - recommended_tag_n / recommended_tag_c
   - purification_methods_ranked  (probabilities or baseline rates)
   - buffer_ph, salt_mm, yield_mg (if trained)
   - suggested_protocol (concise summary)

Python API
----------
```
from predictor import ProteinExpressionPredictor
pred = ProteinExpressionPredictor("./predictor_artifacts")
out = pred.predict("MHHHHHHSSGVDLGTENLYFQSMAS...")
print(out)
```

Notes & Caveats
---------------
- Models use only sequence-derived features (length, composition, hydropathy, crude net charge). 
- Purification method models could not be trained if your dataset only contained positive labels; in that case we rank by empirical success rate.
- Predictions are heuristic and should guide initial screening. Validate experimentally with small-scale expressions.
