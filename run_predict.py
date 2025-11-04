
#!/usr/bin/env python3
import sys, json, os
from predictor import ProteinExpressionPredictor
def main():
    if len(sys.argv) < 2:
        print("Usage: run_predict.py <SEQUENCE or path to FASTA>", file=sys.stderr)
        sys.exit(1)
    arg = sys.argv[1]
    if os.path.exists(arg):
        # Basic FASTA reader: use first sequence
        seq = []
        with open(arg) as f:
            for line in f:
                if line.startswith(">"): 
                    continue
                seq.append(line.strip())
        seq = "".join(seq)
    else:
        seq = arg
    pred = ProteinExpressionPredictor(os.path.dirname(__file__))
    out = pred.predict(seq)
    print(json.dumps(out, indent=2))
if __name__ == "__main__":
    main()
