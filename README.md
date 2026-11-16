# An upper bound for the single-card spectral gap of the overlapping-cycles shuffle

Yitzchak Shmalo

In the overlapping-cycles shuffle of a deck of n cards, the card at position n or the card at position n − k is
moved to the top, each with probability one half. Angel, Peres and Wilson (Ann. Appl. Probab. 18 (2008), 1215–1231)
determined the spectral gap gap(n,k) of the chain followed by a single card, up to a factor 1 + o(1) uniform in k,
as the minimum γ(n,k) of an explicit family of quadratic expressions. They asked whether
gap(n,k) ≤ (2π²/√3 + o(1)) k^(1/2) n^(−2) for all 1 ≤ k < n. The paper shows that it does: γ(n,k) is a constant
multiple of the squared length of the shortest vector of a plane lattice of determinant 2, and Hermite's bound for
plane lattices gives γ(n,k) ≤ (2π²/√3) k^(1/2) n^(−2) for every n. At n = 26, k = 12 the lattice is hexagonal and
the bound is attained.

## Contents

- `paper/main.tex`, `paper/main.pdf` – the paper.
- `code/gap_check.py` – the program of Section 3 (Python 3 with NumPy).
- `SHA256SUMS` – checksums of the files above.

## Reproducing Section 3

```
python code/gap_check.py gamma 2000          # max of gamma(n,k) n^2/sqrt(k) over n <= 2000: 11.396437516 at (26, 12)
python code/gap_check.py eig 300 23          # n = 300, k = 1, 24, 47, ...: max gap n^2/sqrt(k) 8.644; ratio 0.9995 at k = 1
python code/gap_check.py eig 1500 37         # n = 1500, k = 1, 38, 75, ...: max 11.106; ratios 0.791 to 1.025 (k = 260)
python code/gap_check.py conv 0.79 750 1500 3000 6000    # ratios 0.794, 0.791, 0.991, 0.975
python code/gap_check.py conv 0.7407 750 1500 3000       # ratios 0.741, 1.005, 1.002
```

The exact gap is 1 minus the largest modulus of an eigenvalue other than 1 of the n × n transition matrix. The
eigenvalue runs take minutes for n = 1500 and longer for n ≥ 3000.
