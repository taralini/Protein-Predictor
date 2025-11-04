
import streamlit as st
import json, os, io
import pandas as pd
from predictor import ProteinExpressionPredictor

st.set_page_config(page_title="Protein Expression & Purification Predictor", layout="wide")

st.title("🔬 Protein Expression & Purification Predictor")
st.caption("Paste a protein sequence or upload a FASTA. The app predicts expression success, tags, purification methods, and suggests buffer conditions.")

# Sidebar: model/artifact path
default_dir = os.path.join(os.path.dirname(__file__), "predictor_artifacts")
artifact_dir = st.sidebar.text_input("Artifacts directory", value=default_dir)
st.sidebar.info("Artifacts include the trained models and fallback method success rates.")

# Helper: FASTA reader
def read_fasta(file_bytes: bytes) -> str:
    text = file_bytes.decode("utf-8", errors="ignore")
    seq_lines = []
    for line in text.splitlines():
        if line.startswith(">"): 
            # allow multiple entries: join all sequences
            continue
        seq_lines.append(line.strip())
    return "".join(seq_lines)

# Input widgets
tab1, tab2 = st.tabs(["🧬 Paste sequence", "📁 Upload FASTA"])

with tab1:
    seq_text = st.text_area("Protein sequence (one-letter AAs)", height=200, placeholder="MHHHHHHSSGVDLGTENLYFQSMAS...")
with tab2:
    fasta_file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"])
    if fasta_file is not None:
        try:
            seq_text = read_fasta(fasta_file.getvalue())
            st.success(f"Loaded sequence with length {len(seq_text)} from FASTA.")
        except Exception as e:
            st.error(f"Failed to parse FASTA: {e}")

# Run prediction
run = st.button("▶️ Predict", type="primary")
if run:
    try:
        pred = ProteinExpressionPredictor(artifact_dir)
        out = pred.predict(seq_text)

        # Show high-level metrics
        c1, c2, c3 = st.columns(3)
        with c1:
            st.metric("Expression success (prob.)", f"{out.get('expression_success_prob', 0):.0%}")
        with c2:
            st.metric("Buffer pH (suggested)", f"{out.get('buffer_ph', 7.5):.2f}")
        with c3:
            st.metric("Salt (mM, suggested)", f"{out.get('salt_mm', 300):.0f}")

        # Tags
        st.subheader("Recommended tags")
        st.write({
            "N-terminus": out.get("recommended_tag_n"),
            "C-terminus": out.get("recommended_tag_c")
        })

        # Purification methods ranking
        meth = out.get("purification_methods_ranked", [])
        if meth:
            df = pd.DataFrame(meth, columns=["method", "score"])
            df["rank"] = range(1, len(df)+1)
            df = df[["rank", "method", "score"]]
            st.subheader("Purification methods ranked")
            st.dataframe(df, use_container_width=True)
        else:
            st.warning("No per-method model available; consider adding failed attempts (negatives) to your dataset for better discrimination.")

        # Suggested protocol
        st.subheader("Suggested initial protocol")
        st.json(out.get("suggested_protocol", {}))

        # Download JSON
        b = io.BytesIO(json.dumps(out, indent=2).encode("utf-8"))
        st.download_button("Download results (JSON)", data=b, file_name="prediction.json", mime="application/json")

    except Exception as e:
        st.exception(e)

st.markdown("---")
with st.expander("Notes & Caveats"):
    st.markdown("""
- Models use only sequence-derived features (length, composition, hydropathy, crude net charge).
- Purification method rankings may fall back to empirical success rates if your dataset lacks explicit negatives.
- Treat outputs as guidance for first-pass screening; validate with small-scale experiments.
- To improve: add failed attempts as **no** in method columns, expand numeric labels for buffer/yield, and consider host-system metadata.
    """)
