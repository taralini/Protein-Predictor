Streamlit App — Protein Expression & Purification Predictor
===========================================================

Quick start
-----------
1) Create and activate a Python 3.10+ environment.
2) Install deps:
   pip install -r requirements.txt
3) Run:
   streamlit run app.py

Usage
-----
- Paste a protein sequence or upload a FASTA file.
- Click **Predict** to see:
  * Expression success probability
  * Recommended N-/C-terminal tags
  * Ranked purification methods
  * Suggested buffer pH and salt
  * A concise suggested protocol and a JSON export

Customization
-------------
- The app loads models from `./predictor_artifacts` by default.
- Replace that folder with updated artifacts if you retrain on new data.

Improving accuracy
------------------
- Add explicit **no** entries to the `method_success(...)` columns for failed trials.
- Provide more numeric labels for buffer pH, salt, and yields.
