"""Compare a reproduced output folder with the published results in Research_Data/results.
Files must be identical, or equal within a relative difference of 1e-9 (floating-point summation order).
Usage: python3 -I check_against_published.py <reproduced_output_dir> <Research_Data/results>"""
import sys, json, math
from pathlib import Path
import numpy as np, pandas as pd
A, B = Path(sys.argv[1]), Path(sys.argv[2])
same = close = diff = 0; notes = []
def jeq(a, b):
    if isinstance(a, dict): return isinstance(b, dict) and a.keys() == b.keys() and all(jeq(a[k], b[k]) for k in a)
    if isinstance(a, list): return isinstance(b, list) and len(a) == len(b) and all(jeq(p, q) for p, q in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        if a is None or b is None: return a is b
        return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= 1e-9 * max(1, abs(a))
    return a == b
for f in sorted(p for p in B.rglob('*') if p.is_file() and p.suffix in ('.csv', '.json')):
    rel = f.relative_to(B); g = A / rel
    if not g.exists(): notes.append(f'not produced: {rel}'); diff += 1; continue
    if f.read_bytes() == g.read_bytes(): same += 1; continue
    if f.suffix == '.csv':
        x, y = pd.read_csv(f), pd.read_csv(g)
        ok = x.shape == y.shape and list(x.columns) == list(y.columns) and all(
            np.allclose(x[c], y[c], rtol=1e-9, atol=0, equal_nan=True) if pd.api.types.is_numeric_dtype(x[c]) and pd.api.types.is_numeric_dtype(y[c]) else x[c].astype(str).equals(y[c].astype(str)) for c in x.columns)
    else:
        ok = jeq(json.load(open(f)), json.load(open(g)))
    if ok: close += 1; notes.append(f'equal within relative 1e-9 (floating-point order): {rel}')
    else: diff += 1; notes.append(f'DIFFERENT: {rel}')
print(f'identical {same}; equal within 1e-9 {close}; different or missing {diff}')
# analysis datasets supplied in Research_Data/analysis_data, rebuilt by the code in <output>/prepared
AD = B.parent / 'analysis_data'
for f in sorted(AD.glob('*.csv')):
    g = A / 'prepared' / f.name
    print(f'  analysis_data/{f.name}:', 'identical to the rebuilt file' if g.exists() and f.read_bytes() == g.read_bytes() else ('rebuilt file differs' if g.exists() else 'not rebuilt'))
for n in notes: print(' ', n)
