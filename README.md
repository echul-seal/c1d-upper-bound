# C_1d <= 1.455682936505816: certificate and verifiers

Supporting artifact for a new upper bound on C_1d, the signed third autocorrelation constant of the
[Optimization Constants repository](https://teorth.github.io/optimizationproblems/constants/1d.html):
the largest C with max_t |(f*f)(t)| >= C (int f)^2 for every real f supported on [-1/4, 1/4].

## The bound

- Witness: a signed step function with 10240 equal cells on [-1/4, 1/4] and integer cell values
  (largest |v_i| = 10^12, produced by `ub/mkcert.py`).
  - `ub/cons_f10240.json`: key `values`, the input format of the 1d page's `verify.py`.
  - `ub/c_f10240.json`: key `v_int`, the same integers.
- Exact value: max_t |f*f| / (int f)^2 = 76358168487433178963290071040 / 52455219864513498497744894329
  = 1.45568293650581519488... The bound 1.455682936505816 is this rounded up.
- Why a finite computation suffices: with cell width h = 1/20480, f*f is piecewise linear with knots at
  multiples of h, so its max is attained at a knot, where (f*f)(mh) = h * sum_{i+j=m-1} v_i v_j. Hence the ratio is
  2n * max_m |sum_{i+j=m} v_i v_j| / (sum_i v_i)^2 with n = 10240, an exact integer computation.
- Previous merged record: 1.463071582607559* (teorth/optimizationproblems #224).
  Open PRs at the time of writing: #247 (1.460217284022874*), #248 (1.457854884921989*).

## Verify (Python 3, standard library only)

```
sha256sum -c SHA256SUMS
python3 ub/vfy2.py ub/c_f10240.json          # Kronecker substitution, one big-integer product, ~2 s
python3 -I ub/v224.py ub/cons_f10240.json    # the 1d page's verify.py, unchanged (SHA-256 3e75c4ab...0258); exact Fractions, O(n^2), slow
python3 ub/verify_int.py ub/c_f10240.json    # direct O(n^2) integer convolution, slow
```

All three print the fraction above. The `floor25 0.` / `floor 30 digits 0.` lines of `vfy2.py` and `verify_int.py`
print the digits of the ratio (1.4556829...) with a stray `0.` prefix; the fraction line is the result.

## How it was found

L-BFGS on a p-norm smoothing of max|f*f|, p doubled up to 2^23, warm-started by upsampling
1280 -> 2560 -> 5120 -> 10240 cells (`ub/opt1d.py`, logs `ub/l_*.txt`). Resolution trend: 1280 cells 1.456238,
2560 cells 1.455913, 5120 cells 1.455761, 10240 cells 1.455683. The optimized values were scaled and rounded to
integers before certification; the certified number is that of the rounded function.

## Provenance

- Original bundle `c1d_records.tgz`, included byte-for-byte. IPFS CID
  `bafkreifveyzg3tkffwfwwoke4jedupxmkvu2mkr4ep4szg3he3gq7b6ftu`,
  SHA-256 `b526326dcd452d8b6b3944e2483a3eec5569a62a3c23f92c9b6726cd0f87c59d`.
  Gateway: https://gateway.pinata.cloud/ipfs/bafkreifveyzg3tkffwfwwoke4jedupxmkvu2mkr4ep4szg3he3gq7b6ftu
- `NOTES.md` is the bundle's own README.
- `lb/` is a separate, non-record lower-bound experiment from the same bundle and is not part of this claim.

Found and checked with an AI agent (Seal), directed by the human contributor.
