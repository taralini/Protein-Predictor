
import joblib
import numpy as np
import re

KD = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5}
AA = list(KD.keys())
hydrophobic = set(list("AVLIMFWYC"))
aromatic = set(list("FYW"))
polar = set(list("STNQ"))
acidic = set(list("DE"))
basic = set(list("KRH"))

def seq_features(seq: str):
    s = re.sub(r'[^A-Z]', '', str(seq).upper())
    if len(s) == 0:
        raise ValueError("Empty sequence.")
    n = len(s)
    counts = {a: s.count(a) for a in AA}
    gravy = sum(KD.get(a, 0)*counts[a] for a in AA) / n
    feats = {f"aa_{a}_frac": counts[a]/n for a in AA}
    feats.update({
        "length": n,
        "gravy": gravy,
        "frac_hydrophobic": sum(counts.get(a,0) for a in hydrophobic)/n,
        "frac_aromatic": sum(counts.get(a,0) for a in aromatic)/n,
        "frac_polar": sum(counts.get(a,0) for a in polar)/n,
        "frac_acidic": sum(counts.get(a,0) for a in acidic)/n,
        "frac_basic": sum(counts.get(a,0) for a in basic)/n,
        "frac_gly_pro": (counts.get('G',0)+counts.get('P',0))/n,
        "net_charge_pH7": (counts.get('K',0) + counts.get('R',0) + 0.1*counts.get('H',0) - counts.get('D',0) - counts.get('E',0))/n
    })
    keys = ['aa_'+a+'_frac' for a in AA] + ['length','gravy','frac_hydrophobic','frac_aromatic','frac_polar','frac_acidic','frac_basic','frac_gly_pro','net_charge_pH7']
    import numpy as np
    return np.array([feats[k] for k in keys], dtype=float).reshape(1, -1)

class ProteinExpressionPredictor:
    def __init__(self, artifact_dir: str):
        self.expr = joblib.load(artifact_dir + "/clf_expr.pkl")
        self.tag_n = joblib.load(artifact_dir + "/tag_n.pkl")
        self.tag_c = joblib.load(artifact_dir + "/tag_c.pkl")
        self.methods = joblib.load(artifact_dir + "/method_models.pkl")
        import os
        self.reg_ph = joblib.load(artifact_dir + "/reg_buffer_ph.pkl") if os.path.exists(artifact_dir + "/reg_buffer_ph.pkl") else None
        self.reg_salt = joblib.load(artifact_dir + "/reg_salt_mm.pkl") if os.path.exists(artifact_dir + "/reg_salt_mm.pkl") else None
        self.reg_yield = joblib.load(artifact_dir + "/reg_yield_mg.pkl") if os.path.exists(artifact_dir + "/reg_yield_mg.pkl") else None

    def predict(self, sequence: str):
        X = seq_features(sequence)
        out = {}
        out['expression_success_prob'] = float(self.expr.predict_proba(X)[0,1])
        n = self.tag_n.predict(X)[0]
        c = self.tag_c.predict(X)[0]
        out['recommended_tag_n'] = None if n == 'unknown' else str(n)
        out['recommended_tag_c'] = None if c == 'unknown' else str(c)
        m = []
        for name, model in self.methods.items():
            m.append((name, float(model.predict_proba(X)[0,1])))
        m.sort(key=lambda x: x[1], reverse=True)
        out['purification_methods_ranked'] = m
        if self.reg_ph is not None:
            out['buffer_ph'] = float(self.reg_ph.predict(X)[0])
        if self.reg_salt is not None:
            out['salt_mm'] = float(self.reg_salt.predict(X)[0])
        if self.reg_yield is not None:
            out['yield_mg'] = float(self.reg_yield.predict(X)[0])
        # assemble short protocol
        out['suggested_protocol'] = {
            "expression": f"Predicted success: {out['expression_success_prob']:.0%}.",
            "tagging": {"N_term": out['recommended_tag_n'], "C_term": out['recommended_tag_c']},
            "purification": [name for name,_ in m[:3]],
            "buffer": {"pH": out.get('buffer_ph', 7.5), "salt_mM": out.get('salt_mm', 300.0)},
            "expected_yield_mg": out.get('yield_mg', None)
        }
        return out

import json, os
class ProteinExpressionPredictor(ProteinExpressionPredictor):
    def __init__(self, artifact_dir: str):
        super().__init__(artifact_dir)
        self.method_rates = {}
        try:
            with open(artifact_dir + "/method_success_rates.json") as f:
                self.method_rates = json.load(f)
        except Exception:
            self.method_rates = {}

    def predict(self, sequence: str):
        out = super().predict(sequence)
        # If trained models were empty -> backfill from rates
        if len(out['purification_methods_ranked']) == 0 and len(self.method_rates) > 0:
            ranked = sorted(self.method_rates.items(), key=lambda kv: kv[1], reverse=True)
            out['purification_methods_ranked'] = [(k, float(v)) for k,v in ranked]
            out['suggested_protocol']['purification'] = [k for k,_ in ranked[:3]]
        return out
