# C_1d, the signed third autocorrelation constant: new upper and lower bounds

C_1d is the largest C with max_t |(f*f)(t)| >= C (int f)^2 for every real f supported on [-1/4,1/4]
(Optimization constants repository, constants/1d.md). sinc(z) = sin(pi z)/(pi z).

## Upper bound: C_1d <= 1.455682936505816

Resolution trend (same method): 1280 cells 1.456238, 2560 cells 1.455913, 5120 cells 1.455761, 10240 cells 1.455683.

- Witness: signed step function with 10240 equal cells on [-1/4,1/4] and integer cell values:
  ub/cons_f10240.json (key values, the format of the page's verify.py), same values as ub/c_f10240.json (key v_int).
- Exact max|f*f| / (int f)^2 = 76358168487433178963290071040 / 52455219864513498497744894329 = 1.45568293650581519...
  The printed bound is rounded up.
- Previous: 1.463071582607559* (merged, PR 224). Best open PR: 1.457854884921989* (PR 248, 800 cells).
- Three checks, all print the same fraction:
  python -I ub/v224.py ub/cons_f10240.json  (the page's verify.py, unchanged, SHA-256 3e75c4ab...0258)
  python ub/verify_int.py ub/c_f10240.json  (direct integer convolution)
  python ub/vfy2.py ub/c_f10240.json  (Kronecker substitution: one big-integer product)
- Search: L-BFGS on a p-norm smoothing of max|f*f|, p doubled up to 2^23, warm-started by upsampling
  1280 -> 2560 -> 5120 -> 10240 cells (ub/opt1d.py; logs ub/l_*.txt).

## Lower bound (not a record): C_1d >= 1000/831 = 1.2033694344...

The page lists only the trivial bound 1, but github.com/tejasvi/autocorrelation-lower-bounds claims
C_1d >= 400/311 = 1.28617 (Lean-checked), which is stronger. This certificate is an independent, weaker bound by a different method.

Lemma. Scale f so max|f*f| = 1 and suppose (int f)^2 >= 1 - eta. Then e = 1_[-1/2,1/2] - f*f has 0 <= e <= 2 and
int e <= eta, and the Fourier transforms satisfy (f*f)^ = sinc - e^ = (f^)^2, all entire. Take a closed curve G around
z = 1, the only zero of sinc inside it. Suppose at every z = u+iv on G some phase p gives
Re(exp(-ip) sinc z) > sup Re(exp(-ip) e^(z)) over all admissible e (a bathtub: twice the integral of
K(x) = exp(2 pi x v) cos(2 pi x u + p) over the top eta/2 of [-1/2,1/2]). Since t*e is admissible for 0 <= t <= 1,
sinc - t e^ has no zero on G, so (f^)^2 has exactly one zero inside G. Zeros of a square have even multiplicity:
contradiction. Hence (int f)^2 < (1 - eta) max|f*f|, i.e. C_1d >= 1/(1 - eta).

- Certificate: lb/cP_0.169.json, eta = 169/1000: a 180-segment polygon star-shaped about z = 1, from the real axis in
  (1,2) through the upper half plane to the real axis in (0,1). The lower half follows from z -> conj z, p -> -p, e(x) -> e(-x).
- Two independent interval-arithmetic verifiers, both PASS:
  python lb/lbv16.py cP_0.169.json  (python-flint arb, complex balls, 16384 cells; log lb/lbv16.log)
  python lb/verify_lb2.py cP_0.169.json 0 2, then 1 2  (mpmath iv, real-variable formulas, directed rounding,
  4096/32768 cells; logs lb/lb2_0.log, lb/lb2_1.log)
- lb/cert_lb.json is an earlier, weaker certificate (eta 663/5000, C_1d >= 5000/4337 = 1.1529) with its verifiers.

Computed by Seal for Echul, 25 Oct 2026.

## SHA-256
```
c42f5d0733bada1e3d73ed4d704e926b6d6b59d15c1b5589d3c9e829f1b09570  ub/cons_f10240.json
30d501e191e025bc717b47a7a2dd5205aaa35d58eb057e2e551ae9363991a973  ub/c_f10240.json
3e75c4ab0cdb40387765cbd086a572e46605f0c9fd8b8252783ed01aada80258  ub/v224.py
bc7e98e45cd9f75046e95f63af243cfaea54661f8cc4cd43517097f662361fc4  ub/verify_int.py
9f97d899783fe3b54e20b932dc19cb9b1b2e1b4f1334620c146dfb2846581e0f  ub/vfy2.py
723ecd391c25a27e521aafe8ca3e3898833bb636a20c3eca8d9ad5ce6df99d5b  lb/cP_0.169.json
824a74454650ccd0a69aecfca0542f9a3b7f08178c8d7b6d6b71cba3a00f394d  lb/verify_lb2.py
40fd1fc505eb6b2b4b337187a5fb554a78e446aa098b1028c2093dbe1fadde49  lb/lbv16.py
```
