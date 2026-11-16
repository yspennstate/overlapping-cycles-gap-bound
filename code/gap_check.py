"""Numerical checks for the single-card spectral gap of the overlapping-cycles shuffle.

gamma(n,k) = min over m != 0 of (pi^2/(2n^2)) (m^2 k/n + 4n ||m(n-k)/(2n)||^2), with ||x|| = dist(x, Z),
is the quantity of Angel, Peres and Wilson (Ann. Appl. Probab. 18 (2008), formula (1) and Theorem 4).
The bound checked is gamma(n,k) <= (2 pi^2/sqrt 3) sqrt(k)/n^2.

usage: python gap_check.py gamma NMAX       largest gamma(n,k) n^2/sqrt(k) over 2 <= n <= 2000 (and sampled n up to NMAX)
       python gap_check.py eig N STEP       exact gap of the single-card chain against gamma(n,k), k = 1, 1+STEP, ...
       python gap_check.py conv FRAC N...   gap/gamma at k = round(FRAC*N) for each N
"""
import math
import sys

import numpy as np

C = 2 * math.pi ** 2 / math.sqrt(3)


def gamma(n, k):
    A = math.sqrt(k) / n
    mmax = int(math.sqrt(4 / math.sqrt(3) / A)) + 2  # beyond this, A m^2 alone exceeds the bound
    best = float("inf")
    for m in range(1, mmax + 1):
        y = m * (n - k) / (2 * n)
        d = abs(y - round(y))
        best = min(best, (math.pi ** 2 / (2 * n ** 2)) * (m * m * k / n + 4 * n * d * d))
    return best


def kernel(n, k):
    """Transition matrix of the position of one card; position 0 is the top of the deck."""
    P = np.zeros((n, n))
    for i in range(n):
        if i < n - k - 1:
            P[i, i + 1] = 1.0
        elif i == n - k - 1:
            P[i, 0] += 0.5
            P[i, i + 1] += 0.5
        elif i < n - 1:
            P[i, i] += 0.5
            P[i, i + 1] += 0.5
        else:
            P[i, 0] += 0.5
            P[i, i] += 0.5
    return P


def gap(n, k):
    ev = np.linalg.eigvals(kernel(n, k))
    i = int(np.argmin(np.abs(ev - 1.0)))
    return 1.0 - np.delete(np.abs(ev), i).max()


def main():
    mode = sys.argv[1]
    if mode == "gamma":
        nmax = int(sys.argv[2])
        worst = (0.0, None)
        ns = list(range(2, 2001)) + ([int(x) for x in np.geomspace(2001, nmax, 40)] if nmax > 2000 else [])
        for n in ns:
            for k in range(1, n):
                r = gamma(n, k) * n * n / math.sqrt(k)
                if r > worst[0]:
                    worst = (r, (n, k))
        print("max gamma n^2/sqrt(k) = %.9f at %s; bound %.9f; %s"
              % (worst[0], worst[1], C, "ok" if worst[0] <= C * (1 + 1e-12) else "VIOLATION"))
    elif mode == "eig":
        n, step = int(sys.argv[2]), int(sys.argv[3])
        worst = 0.0
        for k in range(1, n, step):
            g, L = gap(n, k), gamma(n, k)
            worst = max(worst, g * n * n / math.sqrt(k))
            print(n, k, "gap %.6e gamma %.6e ratio %.4f" % (g, L, g / L), flush=True)
        print("n", n, "max gap n^2/sqrt(k) %.5f; bound %.5f" % (worst, C))
    elif mode == "conv":
        frac = float(sys.argv[2])
        for n in map(int, sys.argv[3:]):
            k = round(frac * n)
            g, L = gap(n, k), gamma(n, k)
            print(n, k, "gap %.6e gamma %.6e ratio %.4f" % (g, L, g / L), flush=True)


if __name__ == "__main__":
    main()
