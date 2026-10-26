# Optogenetic TdT synthesis simulator v6.1: sensitivity sweep (validated code)

Complete, tested code for a stochastic TdT/light-gated DNA synthesis simulator plus a
sensitivity sweep over **`k_act`, `k_rev`, `dark_leak`, `bias_ratio`**.
Three files, numpy only (pytest for tests; scipy is optional and used by one cross-check).

| file | lines | sha256 (first 12) | purpose |
|---|---:|---|---|
| `tdt_core.py` | 333 | `1d7f7d51daba` | simulator, input validation, metrics, ensemble CLI |
| `tdt_sensitivity.py` | 265 | `7ecd70a20976` | OAT + Latin-hypercube global sweep, statistics, report |
| `test_tdt.py` | 365 | `edd1ee7c3c31` | 54 tests |

## 1. Quick start

```bash
pip install numpy pytest            # scipy optional
# save the three files below next to each other, then:
python -m pytest -v test_tdt.py                                  # ~12 s
python tdt_core.py --target ACGTACGTAC --cells 4 --terminators --set n_strands=300
python tdt_sensitivity.py --mode all --lhs 128 --reps 2 --levels 7 --n-strands 300 \
       --json sweep.json --markdown sweep.md                     # ~1 min on 1 CPU
```

## 2. What was validated

All **54 tests pass** (`54 passed in 11.4s`). The markdown you are reading was also
round-trip checked: the code blocks were extracted from this file into a clean directory and
the suite was re-run there.

| layer | what is checked | against |
|---|---|---|
| input validation | 18 invalid parameter sets are rejected; JSON round-trip of config | explicit expectations |
| helpers | Levenshtein, homopolymer run, truncated product counts as an error | known answers (kitten/sitting = 3, etc.) |
| light model | `activity_step` end/mean activity | RK4 numerical integration of the ODE |
| simulator | Poisson lengths (mean, variance/mean) at constant activity | analytic Poisson limit |
| simulator | mean extension over 3 pulsed cycles | **exact analytic integral of a(t)** (5σ + 0.5 %) |
| simulator | perfect terminators; failed deblock | exact product; Binomial(5, 0.5) deletions |
| simulator | dNTP mass balance, carry-over cannot create molecules, reproducibility | conservation laws, same seed identical |
| statistics | Latin hypercube stratification, log scaling, tie-aware ranks | structural properties, hand-computed Spearman 0.8208, scipy cross-check |
| statistics | SRC, first-order index S1 | synthetic functions with known answers (S1 = 0.8 / 0.2 / 0; non-monotone `(x-0.5)^2`, where Spearman ≈ 0 but S1 > 0.9) |
| physics sanity | more `k_act` / `dark_leak` and lower `k_rev` give more extension; higher `bias_ratio` gives higher fidelity | monotonic directions |
| end to end | full sweeps in each regime, CLI smoke tests | structure and expected dominant parameters |

Two problems the tests caught during development (both fixed, noted for honesty):

1. **The first `free` regime was saturated.** At the inherited `kcat = 8 s^-1` every strand hit the
   capacity cap (`length_ratio` pinned at 3.17), so no parameter could show any effect. The
   regime now uses `kcat = 0.6 s^-1` and a larger capacity (`tests: baseline is not saturated`).
2. **One of my own tests had a wrong expectation** (full-length fraction after failed deblocking):
   the correct law is Binomial(5, 0.5) and the test now checks that, with the simulator unchanged.

## 3. Sweep design

* Parameters and ranges (log-scaled; **placeholders, set from your own priors**):
  `k_act` 0.1-20 s^-1, `k_rev` 0.01-5 s^-1, `dark_leak` 0.001-0.2, `bias_ratio` 1-200.
* **OAT:** 7 levels per parameter, others at baseline (`k_act=2, k_rev=0.5, dark_leak=0.05, bias_ratio=20`); reported as *swing* = max - min of the mean response.
* **Global:** 128-point Latin hypercube, 2 replicate simulations per point, 300 strands, target `ACGTAC`. Per response: Spearman ρ (signed, monotone), SRC (linear), and S1 = share of output variance explained by the parameter alone (8-bin, bias-adjusted, interactions ignored), each with bootstrap 95 % CIs over design points.
* **Responses:** `exact_match` (fraction of strands equal to the target), `edit_per_target_base` (mean Levenshtein distance / target length), `length_ratio` (mean product length / target length).
* **Three regimes**, because the physics differ and one regime hides the others:
  `terminator` (saturated chemistry), `terminator_slow` (`term_kcat_factor = 0.2`, so light dose limits completion), `free` (no terminators, `kcat = 0.6`).
* Replicate noise (mean SD between replicates): `exact_match` ≤ 0.014, `edit_per_target_base` ≤ 0.010, `length_ratio` ≤ 0.014.

## 4. Results (this run, seed 2026)

#### Regime: `terminator`  (target `ACGTAC`, 300 strands, LHS=128, reps=2)

**OAT swing** (max - min of mean response across each parameter's range; baseline k_act=2, k_rev=0.5, dark_leak=0.05, bias_ratio=20)

*exact_match*

| parameter | swing | relative |
|---|---:|---|
| bias_ratio | 0.8700 | `####################` |
| k_act | 0.1117 | `###` |
| k_rev | 0.0600 | `#` |
| dark_leak | 0.0600 | `#` |

*edit_per_target_base*

| parameter | swing | relative |
|---|---:|---|
| bias_ratio | 0.6333 | `####################` |
| k_act | 0.0300 | `#` |
| k_rev | 0.0297 | `#` |
| dark_leak | 0.0186 | `#` |

*length_ratio*

| parameter | swing | relative |
|---|---:|---|
| k_act | 0.0308 | `####################` |
| k_rev | 0.0128 | `########` |
| bias_ratio | 0.0064 | `####` |
| dark_leak | 0.0028 | `##` |

**Global sensitivity** (Latin hypercube; log/linear-scaled unit inputs)

*exact_match*  (linear-model R² = 0.89)

| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |
|---|---|---:|---|
| bias_ratio | +0.97 [+0.94, +0.99] | +0.92 | 0.89 [0.83, 0.96] |
| dark_leak | +0.15 [-0.02, +0.32] | +0.06 | 0.05 [0.00, 0.22] |
| k_act | +0.24 [+0.06, +0.42] | +0.11 | 0.05 [0.00, 0.26] |
| k_rev | +0.01 [-0.17, +0.19] | -0.08 | 0.00 [0.00, 0.17] |

*edit_per_target_base*  (linear-model R² = 0.90)

| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |
|---|---|---:|---|
| bias_ratio | -0.97 [-0.99, -0.94] | -0.94 | 0.92 [0.88, 0.97] |
| dark_leak | -0.14 [-0.31, +0.02] | -0.02 | 0.05 [0.00, 0.21] |
| k_act | -0.24 [-0.42, -0.06] | -0.06 | 0.05 [0.00, 0.26] |
| k_rev | -0.02 [-0.21, +0.15] | +0.05 | 0.03 [0.00, 0.21] |

*length_ratio*  (linear-model R² = 0.29)

| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |
|---|---|---:|---|
| k_rev | -0.48 [-0.61, -0.32] | -0.34 | 0.34 [0.09, 0.49] |
| k_act | +0.14 [-0.05, +0.29] | +0.32 | 0.07 [0.01, 0.26] |
| dark_leak | +0.18 [+0.00, +0.35] | +0.25 | 0.04 [0.00, 0.20] |
| bias_ratio | -0.02 [-0.17, +0.14] | +0.03 | 0.00 [0.00, 0.15] |

#### Regime: `terminator_slow`  (target `ACGTAC`, 300 strands, LHS=128, reps=2)

**OAT swing** (max - min of mean response across each parameter's range; baseline k_act=2, k_rev=0.5, dark_leak=0.05, bias_ratio=20)

*exact_match*

| parameter | swing | relative |
|---|---:|---|
| bias_ratio | 0.5067 | `####################` |
| k_rev | 0.4033 | `################` |
| k_act | 0.3350 | `#############` |
| dark_leak | 0.1433 | `######` |

*edit_per_target_base*

| parameter | swing | relative |
|---|---:|---|
| bias_ratio | 0.5514 | `####################` |
| k_act | 0.3439 | `############` |
| k_rev | 0.3394 | `############` |
| dark_leak | 0.0786 | `###` |

*length_ratio*

| parameter | swing | relative |
|---|---:|---|
| k_act | 0.4356 | `####################` |
| k_rev | 0.4333 | `####################` |
| dark_leak | 0.0983 | `#####` |
| bias_ratio | 0.0442 | `##` |

**Global sensitivity** (Latin hypercube; log/linear-scaled unit inputs)

*exact_match*  (linear-model R² = 0.70)

| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |
|---|---|---:|---|
| bias_ratio | +0.72 [+0.61, +0.81] | +0.66 | 0.43 [0.32, 0.64] |
| k_rev | -0.36 [-0.49, -0.19] | -0.44 | 0.12 [0.07, 0.32] |
| k_act | +0.33 [+0.17, +0.48] | +0.23 | 0.08 [0.03, 0.27] |
| dark_leak | +0.23 [+0.05, +0.40] | +0.11 | 0.06 [0.00, 0.23] |

*edit_per_target_base*  (linear-model R² = 0.74)

| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |
|---|---|---:|---|
| bias_ratio | -0.70 [-0.80, -0.58] | -0.65 | 0.44 [0.33, 0.62] |
| k_rev | +0.36 [+0.19, +0.49] | +0.45 | 0.24 [0.12, 0.43] |
| k_act | -0.35 [-0.50, -0.18] | -0.27 | 0.15 [0.06, 0.33] |
| dark_leak | -0.24 [-0.41, -0.06] | -0.16 | 0.11 [0.00, 0.28] |

*length_ratio*  (linear-model R² = 0.70)

| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |
|---|---|---:|---|
| k_rev | -0.80 [-0.84, -0.74] | -0.68 | 0.61 [0.45, 0.73] |
| k_act | +0.40 [+0.25, +0.51] | +0.41 | 0.11 [0.04, 0.28] |
| dark_leak | +0.20 [+0.03, +0.36] | +0.28 | 0.06 [0.02, 0.27] |
| bias_ratio | -0.06 [-0.23, +0.15] | -0.05 | 0.00 [0.00, 0.14] |

#### Regime: `free`  (target `ACGTAC`, 300 strands, LHS=128, reps=2)

**OAT swing** (max - min of mean response across each parameter's range; baseline k_act=2, k_rev=0.5, dark_leak=0.05, bias_ratio=20)

*exact_match*

| parameter | swing | relative |
|---|---:|---|
| k_act | 0.0033 | `####################` |
| k_rev | 0.0017 | `##########` |
| dark_leak | 0.0017 | `##########` |
| bias_ratio | 0.0017 | `##########` |

*edit_per_target_base*

| parameter | swing | relative |
|---|---:|---|
| k_rev | 1.3231 | `####################` |
| dark_leak | 0.2422 | `####` |
| k_act | 0.1956 | `###` |
| bias_ratio | 0.1939 | `###` |

*length_ratio*

| parameter | swing | relative |
|---|---:|---|
| k_rev | 2.4178 | `####################` |
| k_act | 1.0494 | `#########` |
| dark_leak | 0.6850 | `######` |
| bias_ratio | 0.2108 | `##` |

**Global sensitivity** (Latin hypercube; log/linear-scaled unit inputs)

*exact_match*  (linear-model R² = 0.06)

| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |
|---|---|---:|---|
| k_act | -0.03 [-0.21, +0.16] | -0.04 | 0.03 [0.00, 0.16] |
| bias_ratio | +0.30 [+0.19, +0.43] | +0.25 | 0.03 [0.02, 0.24] |
| k_rev | +0.02 [-0.13, +0.18] | +0.04 | 0.01 [0.00, 0.17] |
| dark_leak | +0.06 [-0.14, +0.23] | -0.03 | 0.00 [0.00, 0.19] |

*edit_per_target_base*  (linear-model R² = 0.67)

| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |
|---|---|---:|---|
| k_rev | -0.60 [-0.70, -0.44] | -0.72 | 0.58 [0.40, 0.71] |
| k_act | +0.29 [+0.09, +0.43] | +0.43 | 0.13 [0.07, 0.33] |
| dark_leak | -0.05 [-0.22, +0.13] | -0.03 | 0.03 [0.00, 0.23] |
| bias_ratio | -0.17 [-0.34, +0.01] | -0.10 | 0.00 [0.00, 0.14] |

*length_ratio*  (linear-model R² = 0.90)

| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |
|---|---|---:|---|
| k_rev | -0.83 [-0.87, -0.77] | -0.84 | 0.70 [0.59, 0.79] |
| k_act | +0.41 [+0.25, +0.52] | +0.47 | 0.16 [0.10, 0.33] |
| dark_leak | +0.21 [+0.05, +0.36] | +0.11 | 0.01 [0.00, 0.21] |
| bias_ratio | -0.06 [-0.23, +0.15] | -0.07 | 0.00 [0.00, 0.16] |


### How to read them

* **Fidelity is set by `bias_ratio` whenever chemistry is saturated.** In `terminator`, `bias_ratio` explains 89-92 % of the variance of `exact_match`/`edit_per_target_base` (ρ = +0.97); the three light parameters sit at the noise floor (S1 ≤ 0.05 with CIs reaching 0). OAT agrees: exact-match swing 0.87 for `bias_ratio` versus ≤ 0.11 for the others.
* **Light parameters matter only when light dose limits the reaction.** In `terminator_slow` they become measurable: `k_rev` (higher reversion, more deletions; S1 0.24 for edit distance), `k_act` (more light, fewer deletions; S1 0.15), `dark_leak` (S1 0.11), while `bias_ratio` is still largest (0.44).
* **Without terminators, extension per cycle is governed by `k_rev`** (S1 = 0.70, ρ = -0.83 for `length_ratio`), then `k_act` (0.16); `bias_ratio` has no effect on length. `exact_match` is essentially 0 everywhere in this regime (largest OAT swing 0.003), so use `edit_per_target_base` and `length_ratio` there. Under these placeholder rates the free regime cannot make correct sequences.
* **`dark_leak` is the weakest of the four in every regime** over 0.001-0.2. A plausible reason is that the slowly decaying tail of light-induced activity (rate `k_rev`) dominates dark-phase activity unless `k_rev` is large, but this mechanism was not tested separately.
* **S1 is a variance share, not an effect size.** E.g. in `terminator`, `k_rev` has S1 = 0.34 for `length_ratio`, yet that response only moves by 0.03 in total (see the swing table). Always read S1 together with swing.

## 5. Limitations (read before trusting numbers)

* **All kinetic defaults are uncalibrated placeholders** (inherited from v5.1 or newly introduced: `kd_strand_uM`, `k_act`, `k_rev`, `dark_leak`, `term_kcat_factor`). The *rankings* above are properties of this model at these ranges, not of any real enzyme/photoreceptor. Re-run with ranges and baselines from your own data.
* Sensitivity depends on the regime: results in one regime do not transfer to another.
* S1 ignores interactions (the linear-model R² of 0.06-0.90 shows how much a linear view captures); a Sobol/Saltelli design would quantify interactions and total-order effects.
* The model is well mixed, one light pulse per cycle at its start, activity held at its exact sub-step mean while events are resolved exactly; strands have no sequence-dependent kinetics or secondary structure.
* Single-core runtime is ~45 s for the full three-regime sweep with the settings above; larger targets or strand counts scale roughly linearly with the number of incorporation events.

## 6. Code

### File: `tdt_core.py`

```python
#!/usr/bin/env python3
"""
tdt_core.py - Optogenetic TdT DNA synthesis simulator (v6.1 core).

MODEL (well-mixed compartment, exact stochastic simulation)
  Per-strand, per-base propensity (s^-1):
      r_i = a(t) * E_bind * kcat_i * s_i / (1 + sum_j s_j) * stall(len) * tf
      s_i = [dNTP_i] / Km_i                       competitive Michaelis-Menten
      E_bind = [TdT] / ([TdT] + Kd_strand)        strand occupancy by enzyme
  Light activation a(t), two-state:
      da/dt = k_act*L(t)*(1-a) - k_rev*(a - dark_leak),   L(t) in {0,1}
  Within a sub-step L is constant, so a(t) is solved exactly and the sub-step
  *mean* of a(t) is used as the rate multiplier (exact in expectation for the
  activity integral). Reaction events inside a sub-step use exact Gillespie.
  Remaining approximation: a(t) is held at its sub-step mean while dNTP
  depletion is resolved event by event (error shrinks with substep_s).

Each cycle: fresh reagent mix for the target base, one light pulse at the start
of the cycle, wash carry-over of leftover dNTP, optional reversible terminator
and deblock.

ALL DEFAULT PARAMETERS ARE PLACEHOLDERS (inherited from v5.1), NOT CALIBRATED.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from dataclasses import asdict, dataclass

import numpy as np

BASES = "ACGT"
MOL_PER_UM_FL = 6.02214076e23 * 1e-15 * 1e-6  # molecules per (uM * fL), ~602
CYCLE_MODES = ("target_only", "biased", "equimolar")


@dataclass
class Params:
    # compartment
    volume_fl: float = 800.0
    n_strands: int = 1000
    # enzyme (s^-1, uM); per-base order A,C,G,T
    tdt_uM: float = 0.8
    kd_strand_uM: float = 0.1
    kcat: tuple = (8.0, 8.0, 8.0, 8.0)
    km_uM: tuple = (100.0, 100.0, 100.0, 100.0)
    # reagents
    dntp_uM: float = 300.0
    cycle_mode: str = "biased"       # target_only | biased | equimolar
    bias_ratio: float = 20.0         # target : each non-target (biased mode)
    wash_carryover: float = 1e-3     # fraction of leftover dNTP carried over
    # light / activation
    cycle_s: float = 10.0
    pulse_on_s: float = 1.0          # one pulse at the start of each cycle
    substep_s: float = 0.25          # cycle_s and pulse_on_s must be multiples
    k_act: float = 2.0               # s^-1 under light
    k_rev: float = 0.5               # s^-1 dark reversion
    dark_leak: float = 0.05          # residual active fraction in the dark
    # chemistry
    terminators: bool = False
    term_eff: float = 0.99
    deblock_eff: float = 1.0
    term_kcat_factor: float = 1.0    # terminators are usually much slower
    # stalling of long strands
    stall_len: int = 15
    stall_factor: float = 0.2
    max_len_factor: float = 3.0      # strand capacity = factor * target length

    def validate(self) -> "Params":
        e = []

        def req(cond, msg):
            if not cond:
                e.append(msg)

        req(self.volume_fl > 0, "volume_fl must be > 0")
        req(int(self.n_strands) >= 1, "n_strands must be >= 1")
        req(self.tdt_uM > 0, "tdt_uM must be > 0")
        req(self.kd_strand_uM >= 0, "kd_strand_uM must be >= 0")
        req(len(self.kcat) == 4 and all(k >= 0 for k in self.kcat), "kcat needs 4 values >= 0")
        req(len(self.km_uM) == 4 and all(k > 0 for k in self.km_uM), "km_uM needs 4 values > 0")
        req(self.dntp_uM > 0, "dntp_uM must be > 0")
        req(self.cycle_mode in CYCLE_MODES, f"cycle_mode must be one of {CYCLE_MODES}")
        req(self.bias_ratio >= 1, "bias_ratio must be >= 1")
        req(0 <= self.wash_carryover <= 1, "wash_carryover must be in [0,1]")
        req(self.cycle_s > 0 and self.substep_s > 0, "cycle_s and substep_s must be > 0")
        req(0 <= self.pulse_on_s <= self.cycle_s, "pulse_on_s must be in [0, cycle_s]")
        if self.cycle_s > 0 and self.substep_s > 0:
            for nm, v in (("cycle_s", self.cycle_s), ("pulse_on_s", self.pulse_on_s)):
                q = v / self.substep_s
                req(abs(q - round(q)) < 1e-6, f"{nm} must be a multiple of substep_s")
        req(self.k_act >= 0 and self.k_rev >= 0, "k_act and k_rev must be >= 0")
        req(0 <= self.dark_leak <= 1, "dark_leak must be in [0,1]")
        req(0 <= self.term_eff <= 1 and 0 <= self.deblock_eff <= 1,
            "term_eff and deblock_eff must be in [0,1]")
        req(self.term_kcat_factor >= 0, "term_kcat_factor must be >= 0")
        req(int(self.stall_len) >= 1, "stall_len must be >= 1")
        req(0 <= self.stall_factor <= 1, "stall_factor must be in [0,1]")
        req(self.max_len_factor >= 1, "max_len_factor must be >= 1")
        if e:
            raise ValueError("; ".join(e))
        return self


def _round(x, rng):
    """Stochastic rounding: keeps the expected molecule count unbiased."""
    fl = np.floor(x)
    return fl + (rng.random(x.shape) < (x - fl))


def activity_step(a, light, p: Params):
    """Exact two-state solution over one sub-step.
    Returns (a_end, mean_a_over_substep)."""
    kl = p.k_act if light else 0.0
    den = kl + p.k_rev
    if den == 0.0:
        return a, a
    a_inf = (kl + p.k_rev * p.dark_leak) / den
    x = den * p.substep_s
    decay = np.exp(-x)
    a_end = a_inf + (a - a_inf) * decay
    mean = a_inf + (a - a_inf) * (1.0 - decay) / x
    return a_end, mean


def simulate_cell(p: Params, target: str, seed) -> dict:
    p.validate()
    if not target or any(c not in BASES for c in target):
        raise ValueError("target must be a non-empty string over ACGT")
    rng = np.random.default_rng(seed)
    L, N = len(target), int(p.n_strands)
    max_len = int(L * p.max_len_factor) + 1
    seqs = np.full((N, max_len), -1, np.int8)
    length = np.zeros(N, np.int64)
    blocked = np.zeros(N, bool)
    kcat, km = np.asarray(p.kcat, float), np.asarray(p.km_uM, float)
    mol_per_uM = p.volume_fl * MOL_PER_UM_FL
    e_bind = p.tdt_uM / (p.tdt_uM + p.kd_strand_uM)
    tf = p.term_kcat_factor if p.terminators else 1.0
    dt = p.substep_s
    n_sub = int(round(p.cycle_s / dt))
    a = p.dark_leak
    carry = np.zeros(4)
    supplied = consumed = 0

    for k in range(L):
        t_idx = BASES.index(target[k])
        conc = np.zeros(4)
        if p.cycle_mode == "equimolar":
            conc[:] = p.dntp_uM
        elif p.cycle_mode == "biased":
            conc[:] = p.dntp_uM / p.bias_ratio
            conc[t_idx] = p.dntp_uM
        else:  # target_only
            conc[t_idx] = p.dntp_uM
        fresh = _round(conc * mol_per_uM, rng)
        supplied += int(fresh.sum())
        mol = fresh + carry

        if p.terminators:
            blocked[blocked & (rng.random(N) < p.deblock_eff)] = False
        w = np.where(blocked | (length >= max_len), 0.0,
                     np.where(length >= p.stall_len, p.stall_factor, 1.0))

        for s in range(n_sub):
            light = (s * dt) < p.pulse_on_s - 1e-9
            a, f = activity_step(a, light, p)
            tl = 0.0
            while True:
                sat = (mol / mol_per_uM) / km
                rb = f * e_bind * tf * kcat * sat / (1.0 + sat.sum())
                r1, wsum = rb.sum(), w.sum()
                if r1 <= 0 or wsum <= 0:
                    break
                tl += rng.exponential(1.0 / (r1 * wsum))
                if tl > dt:
                    break
                cs = np.cumsum(w)
                i = int(np.searchsorted(cs, rng.random() * cs[-1], side="right"))
                b = min(int(np.searchsorted(np.cumsum(rb), rng.random() * r1, side="right")), 3)
                seqs[i, length[i]] = b
                length[i] += 1
                mol[b] -= 1
                consumed += 1
                if p.terminators and rng.random() < p.term_eff:
                    blocked[i] = True
                w[i] = 0.0 if (blocked[i] or length[i] >= max_len) else (
                    p.stall_factor if length[i] >= p.stall_len else 1.0)

        assert mol.min() >= 0, "negative dNTP count"
        carry = _round(mol * p.wash_carryover, rng)

    return dict(seqs=seqs, length=length, supplied=supplied, consumed=consumed)


# ---------------------------- analysis ----------------------------------------
def levenshtein(a, b) -> int:
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


def max_run(row) -> int:
    best = cur = 1 if len(row) else 0
    for i in range(1, len(row)):
        cur = cur + 1 if row[i] == row[i - 1] else 1
        best = max(best, cur)
    return best


def analyse(res: dict, target: str) -> dict:
    tgt = [BASES.index(c) for c in target]
    L = len(tgt)
    seqs, length = res["seqs"], res["length"]
    N = len(length)
    edits = np.empty(N)
    runs = np.empty(N)
    for n in range(N):
        row = seqs[n, :length[n]].tolist()
        edits[n] = 0 if row == tgt else levenshtein(row, tgt)
        runs[n] = max_run(row)
    full = length == L
    pos_err = ((seqs[full, :L] != np.array(tgt)).mean(axis=0).tolist()
               if full.any() else [float("nan")] * L)
    return dict(
        exact_match=float((edits == 0).mean()),
        full_length=float(full.mean()),
        mean_edit_distance=float(edits.mean()),
        edit_per_target_base=float(edits.mean() / L),
        mean_length=float(length.mean()),
        length_ratio=float(length.mean() / L),
        mean_max_homopolymer=float(runs.mean()),
        mean_dntp_use_per_strand=float(res["consumed"] / N),
        per_position_mismatch_full_length=pos_err,
    )


# ---------------------------- ensemble ----------------------------------------
def _worker(args):
    pdict, target, seed, fasta_n = args
    p = Params(**pdict)
    res = simulate_cell(p, target, seed)
    out = analyse(res, target)
    if fasta_n:
        out["_fasta"] = ["".join(BASES[x] for x in res["seqs"][n, :res["length"][n]])
                         for n in range(min(fasta_n, int(p.n_strands)))]
    return out


def bootstrap_ci(x, rng, n=2000):
    x = np.asarray(x, float)
    means = rng.choice(x, size=(n, len(x)), replace=True).mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def run_ensemble(p: Params, target: str, cells: int, seed: int, procs: int, fasta_n=0):
    seeds = [int(s.generate_state(1)[0]) for s in np.random.SeedSequence(seed).spawn(cells)]
    jobs = [(asdict(p), target, s, fasta_n if i == 0 else 0) for i, s in enumerate(seeds)]
    if procs > 1 and cells > 1:
        with mp.Pool(procs) as pool:
            outs = pool.map(_worker, jobs)
    else:
        outs = [_worker(j) for j in jobs]
    rng = np.random.default_rng(seed)
    summary = {}
    for key in outs[0]:
        if key in ("per_position_mismatch_full_length", "_fasta"):
            continue
        v = [o[key] for o in outs]
        lo, hi = bootstrap_ci(v, rng) if cells > 1 else (v[0], v[0])
        summary[key] = dict(mean=float(np.mean(v)), ci95=[lo, hi])
    summary["per_position_mismatch_full_length"] = np.nanmean(
        [o["per_position_mismatch_full_length"] for o in outs], axis=0).tolist()
    return summary, outs[0].get("_fasta", [])


def parse_set(items):
    out = {}
    for it in items or []:
        k, v = it.split("=", 1)
        try:
            v = json.loads(v)
        except json.JSONDecodeError:
            pass
        out[k] = tuple(v) if isinstance(v, list) else v
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Optogenetic TdT synthesis simulator")
    ap.add_argument("--target", default="ACGTACGTAC")
    ap.add_argument("--cells", type=int, default=4)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--procs", type=int, default=max(1, mp.cpu_count() - 1))
    ap.add_argument("--terminators", action="store_true")
    ap.add_argument("--set", action="append", metavar="KEY=VALUE", help="override any Params field")
    ap.add_argument("--out", help="write config + results JSON here")
    ap.add_argument("--fasta", help="write first cell's products to FASTA")
    ap.add_argument("--fasta-n", type=int, default=100)
    a = ap.parse_args(argv)

    overrides = parse_set(a.set)
    unknown = set(overrides) - set(Params.__dataclass_fields__)
    if unknown:
        raise SystemExit(f"unknown parameter(s): {sorted(unknown)}")
    p = Params(**overrides)
    if a.terminators:
        p.terminators = True
    p.validate()
    summary, fasta = run_ensemble(p, a.target.upper(), a.cells, a.seed, a.procs,
                                  a.fasta_n if a.fasta else 0)
    for k, v in summary.items():
        if isinstance(v, dict):
            print(f"{k:32s} {v['mean']:.4g}  95% CI [{v['ci95'][0]:.4g}, {v['ci95'][1]:.4g}]")
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(dict(params=asdict(p), target=a.target, cells=a.cells, seed=a.seed,
                           numpy=np.__version__, results=summary), fh, indent=2)
    if a.fasta:
        with open(a.fasta, "w") as fh:
            for i, s in enumerate(fasta):
                fh.write(f">strand_{i}\n{s}\n")


if __name__ == "__main__":
    main()
```

### File: `tdt_sensitivity.py`

```python
#!/usr/bin/env python3
"""
tdt_sensitivity.py - Sensitivity sweep over k_act, k_rev, dark_leak, bias_ratio.

Two complementary analyses (numpy only, no SALib/scipy needed):
  1. OAT  (one-at-a-time): each parameter swept over its range with the others
     at baseline -> response curves and "swing" (max - min of the mean response).
  2. GLOBAL: Latin-hypercube design over the whole space, then per response
        - Spearman rank correlation (monotone effect, signed)
        - SRC: standardized regression coefficients on the unit (log/linear
          scaled) inputs, with the linear-model R^2
        - S1: first-order variance share Var(E[Y|X_i])/Var(Y) from quantile
          bins, adjusted for the bins' null bias (detects NON-monotone effects
          that Spearman misses). Approximate; it ignores interactions.
     with bootstrap 95% CIs over design points.

Ranges in PARAM_SPACE are PLACEHOLDERS: set them from your own prior knowledge.
Three regimes are swept because the physics differ:
  terminator      : saturated chemistry; light parameters only change speed
  terminator_slow : slow terminator incorporation; light dose limits completion
  free            : no terminators; light dose controls how many bases are added
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from dataclasses import asdict

import numpy as np

from tdt_core import Params, _worker

PARAM_SPACE = {            # name: (low, high, scale)
    "k_act": (0.1, 20.0, "log"),
    "k_rev": (0.01, 5.0, "log"),
    "dark_leak": (0.001, 0.2, "log"),
    "bias_ratio": (1.0, 200.0, "log"),
}
NAMES = tuple(PARAM_SPACE)
RESPONSES = ("exact_match", "edit_per_target_base", "length_ratio")
# Regimes are tuned so the BASELINE is not saturated; otherwise no parameter can
# show an effect (at kcat=8 s^-1 a 10 s cycle adds ~12 bases/strand and every strand
# hits the capacity cap).  kcat / term_kcat_factor are placeholders: calibrate them.
MODES = {
    # saturated one-base-per-cycle chemistry: only reagent selectivity should matter
    "terminator": dict(terminators=True, term_eff=0.99),
    # slow terminator incorporation: light dose limits per-cycle completion (deletions)
    "terminator_slow": dict(terminators=True, term_eff=0.99, term_kcat_factor=0.2),
    # no terminators: light dose sets how many bases are added (~1/cycle at baseline)
    "free": dict(terminators=False, kcat=(0.6, 0.6, 0.6, 0.6), max_len_factor=8.0),
}


# ------------------------------ design ----------------------------------------
def latin_hypercube(n: int, d: int, rng) -> np.ndarray:
    """n points in [0,1)^d; each column has exactly one point per 1/n stratum."""
    u = np.empty((n, d))
    for j in range(d):
        u[:, j] = (rng.permutation(n) + rng.random(n)) / n
    return u


def scale_unit(u: np.ndarray, names=NAMES, space=PARAM_SPACE) -> np.ndarray:
    x = np.empty_like(u, dtype=float)
    for j, nm in enumerate(names):
        lo, hi, sc = space[nm]
        x[:, j] = lo * (hi / lo) ** u[:, j] if sc == "log" else lo + (hi - lo) * u[:, j]
    return x


def evaluate(points, names, base: Params, target, reps, seed, procs):
    """Run the simulator at each row of `points`. Returns (mean, sd) arrays of
    shape (n_points, len(RESPONSES)), averaged over `reps` replicate runs."""
    n = len(points)
    seeds = np.random.SeedSequence(seed).spawn(n * reps)
    jobs = []
    for i, row in enumerate(points):
        pd = asdict(base)
        pd.update({nm: float(v) for nm, v in zip(names, row)})
        for r in range(reps):
            jobs.append((pd, target, int(seeds[i * reps + r].generate_state(1)[0]), 0))
    if procs > 1 and len(jobs) > 1:
        with mp.Pool(procs) as pool:
            outs = pool.map(_worker, jobs)
    else:
        outs = [_worker(j) for j in jobs]
    Y = np.array([[o[k] for k in RESPONSES] for o in outs]).reshape(n, reps, len(RESPONSES))
    sd = Y.std(axis=1, ddof=1) if reps > 1 else np.zeros((n, len(RESPONSES)))
    return Y.mean(axis=1), sd


# ------------------------------ statistics ------------------------------------
def rank(x) -> np.ndarray:
    """1-based ranks with average ranks for ties."""
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    ranks = np.empty(len(x))
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and xs[j + 1] == xs[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def spearman(x, y) -> float:
    rx, ry = rank(x), rank(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def src(X, y):
    """Standardized regression coefficients and R^2 of the linear fit."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    sx, sy = X.std(axis=0), y.std()
    if sy == 0 or np.any(sx == 0):
        return np.full(X.shape[1], np.nan), float("nan")
    Xs, ys = (X - X.mean(axis=0)) / sx, (y - y.mean()) / sy
    beta, *_ = np.linalg.lstsq(Xs, ys, rcond=None)
    r2 = 1.0 - np.sum((ys - Xs @ beta) ** 2) / np.sum(ys ** 2)
    return beta, float(r2)


def first_order_index(x, y, bins: int = 8) -> float:
    """Adjusted between-bin variance share of y (quantile bins of x), clipped at 0."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(y)
    vy = y.var()
    if vy == 0 or n <= bins:
        return float("nan")
    edges = np.quantile(x, np.linspace(0, 1, bins + 1)[1:-1])
    idx = np.searchsorted(edges, x, side="right")
    ybar = y.mean()
    between = sum(np.sum(idx == b) * (y[idx == b].mean() - ybar) ** 2
                  for b in range(bins) if np.any(idx == b)) / n
    s = between / vy
    adj = 1.0 - (1.0 - s) * (n - 1) / (n - bins)
    return float(max(adj, 0.0))


def global_indices(U, Y, names=NAMES, responses=RESPONSES, rng=None, n_boot=300, bins=8):
    """Per response and parameter: Spearman, SRC, S1 (+ bootstrap CIs)."""
    rng = rng or np.random.default_rng(0)
    n = len(U)
    out = {}
    boot = [rng.integers(0, n, n) for _ in range(n_boot)]
    for k, resp in enumerate(responses):
        y = Y[:, k]
        beta, r2 = src(U, y)
        entry = {"linear_R2": r2, "params": {}}
        for j, nm in enumerate(names):
            rho = spearman(U[:, j], y)
            s1 = first_order_index(U[:, j], y, bins)
            rb = [spearman(U[b, j], y[b]) for b in boot]
            sb = [first_order_index(U[b, j], y[b], bins) for b in boot]
            entry["params"][nm] = dict(
                spearman=rho, spearman_ci95=[float(np.nanpercentile(rb, 2.5)), float(np.nanpercentile(rb, 97.5))],
                src=float(beta[j]), s1=s1,
                s1_ci95=[float(np.nanpercentile(sb, 2.5)), float(np.nanpercentile(sb, 97.5))])
        out[resp] = entry
    return out


# ------------------------------ sweeps ----------------------------------------
def oat_sweep(base: Params, target, levels, reps, seed, procs, names=NAMES, space=PARAM_SPACE):
    baseline = np.array([getattr(base, nm) for nm in names], float)
    result = {}
    for j, nm in enumerate(names):
        lo, hi, sc = space[nm]
        grid = np.geomspace(lo, hi, levels) if sc == "log" else np.linspace(lo, hi, levels)
        pts = np.tile(baseline, (levels, 1))
        pts[:, j] = grid
        m, sd = evaluate(pts, names, base, target, reps, seed + 1000 * (j + 1), procs)
        result[nm] = dict(
            grid=grid.tolist(),
            mean={r: m[:, k].tolist() for k, r in enumerate(RESPONSES)},
            sd={r: sd[:, k].tolist() for k, r in enumerate(RESPONSES)},
            swing={r: float(m[:, k].max() - m[:, k].min()) for k, r in enumerate(RESPONSES)})
    return dict(baseline=dict(zip(names, baseline.tolist())), params=result)


def run_sweep(mode, target="ACGTAC", n_strands=300, lhs=128, reps=2, levels=7,
              seed=2026, procs=1, n_boot=300, base_overrides=None):
    if mode not in MODES:
        raise ValueError(f"mode must be one of {tuple(MODES)}")
    base = Params(n_strands=n_strands, **MODES[mode], **(base_overrides or {})).validate()
    rng = np.random.default_rng(seed)
    U = latin_hypercube(lhs, len(NAMES), rng)
    Y, Ysd = evaluate(scale_unit(U), NAMES, base, target, reps, seed, procs)
    return dict(
        mode=mode, target=target, base=asdict(base),
        design=dict(lhs=lhs, reps=reps, levels=levels, seed=seed, n_boot=n_boot, space=PARAM_SPACE),
        oat=oat_sweep(base, target, levels, reps, seed, procs),
        global_=global_indices(U, Y, rng=rng, n_boot=n_boot),
        mean_replicate_sd={r: float(Ysd[:, k].mean()) for k, r in enumerate(RESPONSES)})


# ------------------------------ reporting -------------------------------------
def _bar(v, width=20):
    v = 0.0 if v != v else max(0.0, min(1.0, v))
    return "#" * int(round(v * width))


def to_markdown(res: dict) -> str:
    L = [f"### Regime: `{res['mode']}`  (target `{res['target']}`, "
         f"{res['base']['n_strands']} strands, LHS={res['design']['lhs']}, "
         f"reps={res['design']['reps']})", ""]
    L += ["**OAT swing** (max - min of mean response across each parameter's range; "
          "baseline " + ", ".join(f"{k}={v:g}" for k, v in res["oat"]["baseline"].items()) + ")", ""]
    for resp in RESPONSES:
        sw = {nm: res["oat"]["params"][nm]["swing"][resp] for nm in NAMES}
        top = max(sw.values()) or 1.0
        L += [f"*{resp}*", "", "| parameter | swing | relative |", "|---|---:|---|"]
        for nm in sorted(sw, key=sw.get, reverse=True):
            L.append(f"| {nm} | {sw[nm]:.4f} | `{_bar(sw[nm] / top)}` |")
        L.append("")
    L += ["**Global sensitivity** (Latin hypercube; log/linear-scaled unit inputs)", ""]
    for resp in RESPONSES:
        g = res["global_"][resp]
        L += [f"*{resp}*  (linear-model R² = {g['linear_R2']:.2f})", "",
              "| parameter | Spearman ρ [95% CI] | SRC | S1 [95% CI] |", "|---|---|---:|---|"]
        order = sorted(NAMES, key=lambda n: -(g["params"][n]["s1"] if g["params"][n]["s1"] == g["params"][n]["s1"] else 0))
        for nm in order:
            q = g["params"][nm]
            L.append(f"| {nm} | {q['spearman']:+.2f} [{q['spearman_ci95'][0]:+.2f}, {q['spearman_ci95'][1]:+.2f}] "
                     f"| {q['src']:+.2f} | {q['s1']:.2f} [{q['s1_ci95'][0]:.2f}, {q['s1_ci95'][1]:.2f}] |")
        L.append("")
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Sensitivity sweep for the TdT simulator")
    ap.add_argument("--mode", choices=list(MODES) + ["all"], default="all")
    ap.add_argument("--target", default="ACGTAC")
    ap.add_argument("--n-strands", type=int, default=300)
    ap.add_argument("--lhs", type=int, default=128, help="Latin-hypercube design points")
    ap.add_argument("--reps", type=int, default=2, help="replicate simulations per point")
    ap.add_argument("--levels", type=int, default=7, help="OAT grid levels per parameter")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--procs", type=int, default=max(1, mp.cpu_count() - 1))
    ap.add_argument("--json", help="write full results JSON")
    ap.add_argument("--markdown", help="write report tables (markdown)")
    a = ap.parse_args(argv)

    modes = list(MODES) if a.mode == "all" else [a.mode]
    results = [run_sweep(m, a.target.upper(), a.n_strands, a.lhs, a.reps, a.levels, a.seed, a.procs)
               for m in modes]
    md = "\n".join(to_markdown(r) for r in results)
    print(md)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(results, fh, indent=2)
    if a.markdown:
        with open(a.markdown, "w") as fh:
            fh.write(md)


if __name__ == "__main__":
    main()
```

### File: `test_tdt.py`

```python
"""
Test suite for tdt_core.py and tdt_sensitivity.py.   Run:  python -m pytest -v test_tdt.py

Layers:
  1. input validation and helpers (known-answer tests)
  2. simulator vs analytic limits (Poisson, activity integral, terminators, mass balance)
  3. statistics code vs hand/known answers and synthetic functions with known indices
  4. physical sanity (monotonic responses) and end-to-end sweep
All stochastic tests use fixed seeds, so they are deterministic.
"""
import json
from dataclasses import replace

import numpy as np
import pytest

import tdt_core as core
import tdt_sensitivity as sens
from tdt_core import Params, activity_step, analyse, levenshtein, max_run, simulate_cell


# ============================ 1. validation & helpers ===========================
@pytest.mark.parametrize("bad", [
    dict(volume_fl=0), dict(n_strands=0), dict(tdt_uM=0), dict(kcat=(1, 2, 3)),
    dict(km_uM=(1, 1, 1, 0)), dict(dntp_uM=-1), dict(cycle_mode="nope"),
    dict(bias_ratio=0.5), dict(wash_carryover=1.5), dict(dark_leak=1.2),
    dict(k_act=-1), dict(k_rev=-0.1), dict(term_eff=2), dict(stall_factor=-0.1),
    dict(pulse_on_s=20.0), dict(cycle_s=10.1), dict(substep_s=0.3), dict(max_len_factor=0.5),
])
def test_params_validation_rejects(bad):
    with pytest.raises(ValueError):
        Params(**bad).validate()


def test_params_defaults_valid_and_json_roundtrip():
    p = Params().validate()
    from dataclasses import asdict
    q = Params(**json.loads(json.dumps(asdict(p))))
    assert q.validate().k_act == p.k_act


def test_bad_target_rejected():
    for t in ("", "ACGX", "acg"):
        with pytest.raises(ValueError):
            simulate_cell(Params(n_strands=5), t, 0)


def test_levenshtein_known_answers():
    assert levenshtein(list("kitten"), list("sitting")) == 3
    assert levenshtein([], [1, 2, 3]) == 3
    assert levenshtein([1, 2, 3], [1, 2, 3]) == 0
    assert levenshtein([0, 1, 2, 3], [0, 2, 3]) == 1          # one deletion
    assert levenshtein([0, 1], [0, 1, 1]) == 1                # one insertion
    assert levenshtein([0, 1, 2], [0, 3, 2]) == 1             # one substitution


def test_truncated_product_counts_as_error():
    """The old Hamming/zip metric scored truncations as perfect; edit distance must not."""
    res = dict(seqs=np.array([[0, 1, -1, -1]], np.int8), length=np.array([2]), consumed=2)
    out = analyse(res, "ACGT")
    assert out["exact_match"] == 0.0 and out["mean_edit_distance"] == 2.0


def test_max_run():
    assert max_run([]) == 0 and max_run([1]) == 1
    assert max_run([0, 0, 1, 1, 1, 0]) == 3


def test_activity_step_matches_numeric_ode():
    p = Params(k_act=3.0, k_rev=0.7, dark_leak=0.1, substep_s=0.25)
    for light in (True, False):
        a0 = 0.3
        a_end, mean = activity_step(a0, light, p)
        kl = p.k_act if light else 0.0
        a, tot, n = a0, 0.0, 20000
        h = p.substep_s / n

        def f(a_):
            return kl * (1 - a_) - p.k_rev * (a_ - p.dark_leak)
        for _ in range(n):                       # RK4
            k1 = f(a); k2 = f(a + h * k1 / 2); k3 = f(a + h * k2 / 2); k4 = f(a + h * k3)
            a_new = a + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            tot += 0.5 * (a + a_new) * h
            a = a_new
        assert a_end == pytest.approx(a, abs=1e-8)
        assert mean == pytest.approx(tot / p.substep_s, abs=1e-6)


def test_activity_step_zero_rates_is_identity():
    p = Params(k_act=0.0, k_rev=0.0)
    assert activity_step(0.4, False, p) == (0.4, 0.4)


# ============================ 2. simulator vs analytics =========================
def analytic_activity_integral(p: Params, n_cycles: int) -> float:
    """Exact integral of a(t) over n cycles (light for pulse_on_s, dark afterwards)."""
    a, total = p.dark_leak, 0.0
    for _ in range(n_cycles):
        for kl, T in ((p.k_act, p.pulse_on_s), (0.0, p.cycle_s - p.pulse_on_s)):
            den = kl + p.k_rev
            if T == 0:
                continue
            if den == 0:
                total += a * T
                continue
            a_inf = (kl + p.k_rev * p.dark_leak) / den
            total += a_inf * T + (a - a_inf) * (1 - np.exp(-den * T)) / den
            a = a_inf + (a - a_inf) * np.exp(-den * T)
    return total


NO_DEPLETION = dict(volume_fl=1e6, cycle_mode="target_only", stall_len=10 ** 9,
                    max_len_factor=100, wash_carryover=0.0)


def test_poisson_limit_constant_activity():
    p = Params(n_strands=500, k_act=0.0, k_rev=1.0, dark_leak=1.0, **NO_DEPLETION)
    res = simulate_cell(p, "A", seed=1)
    sat = p.dntp_uM / p.km_uM[0]
    lam = p.tdt_uM / (p.tdt_uM + p.kd_strand_uM) * p.kcat[0] * sat / (1 + sat) * p.cycle_s
    m, v = res["length"].mean(), res["length"].var()
    assert abs(m - lam) < 5 * np.sqrt(lam / 500)
    assert abs(v / lam - 1) < 0.25


def test_mean_extension_matches_analytic_activity_integral():
    """Validates the light model + sub-step scheme end to end (3 pulsed cycles)."""
    p = Params(n_strands=800, k_act=2.5, k_rev=0.6, dark_leak=0.03, **NO_DEPLETION)
    res = simulate_cell(p, "AAA", seed=11)
    sat = p.dntp_uM / p.km_uM[0]
    rate = p.tdt_uM / (p.tdt_uM + p.kd_strand_uM) * p.kcat[0] * sat / (1 + sat)
    expected = rate * analytic_activity_integral(p, 3)
    got = res["length"].mean()
    assert abs(got - expected) < 5 * np.sqrt(expected / 800) + 0.005 * expected, (got, expected)


def test_perfect_terminators_give_exact_product():
    p = Params(n_strands=200, volume_fl=1e5, cycle_mode="target_only", terminators=True,
               term_eff=1.0, deblock_eff=1.0, wash_carryover=0.0)
    out = analyse(simulate_cell(p, "ACGTTA", seed=2), "ACGTTA")
    assert out["exact_match"] == 1.0 and out["mean_max_homopolymer"] >= 2


def test_failed_deblock_creates_deletions_not_insertions():
    p = Params(n_strands=300, volume_fl=1e5, cycle_mode="target_only", terminators=True,
               term_eff=1.0, deblock_eff=0.5, wash_carryover=0.0)
    r = simulate_cell(p, "ACGTAC", seed=4)
    L = r["length"]
    assert L.max() <= 6 and L.min() >= 1           # never longer than the target
    # Cycle 1 always extends (no blocked strands yet). In each of the 5 later cycles every
    # strand is blocked and unblocks independently with p=0.5, so L-1 ~ Binomial(5, 0.5).
    se = np.sqrt(5 * 0.25 / 300)
    assert abs((L - 1).mean() - 2.5) < 5 * se
    assert abs((L == 6).mean() - 0.5 ** 5) < 5 * np.sqrt(0.5 ** 5 * (1 - 0.5 ** 5) / 300)


def test_mass_conservation_under_depletion():
    p = Params(n_strands=500, volume_fl=0.002, cycle_mode="target_only", wash_carryover=0.0)
    r = simulate_cell(p, "AC", seed=3)
    assert r["consumed"] == r["length"].sum() <= r["supplied"]
    assert r["consumed"] >= 0.9 * r["supplied"]     # dNTP is the limiting reagent


def test_carryover_never_creates_molecules():
    p = Params(n_strands=200, volume_fl=0.01, cycle_mode="biased", wash_carryover=1.0)
    r = simulate_cell(p, "ACGT", seed=5)
    assert r["consumed"] <= r["supplied"]


def test_reproducible_and_seed_sensitive():
    p = Params(n_strands=80, terminators=True, volume_fl=1e4)
    a = simulate_cell(p, "ACGTAC", 7)["seqs"]
    b = simulate_cell(p, "ACGTAC", 7)["seqs"]
    c = simulate_cell(p, "ACGTAC", 8)["seqs"]
    assert np.array_equal(a, b) and not np.array_equal(a, c)


def test_analyse_outputs_in_range():
    p = Params(n_strands=100, terminators=True, volume_fl=1e4)
    out = analyse(simulate_cell(p, "ACGTAC", 9), "ACGTAC")
    for k in ("exact_match", "full_length"):
        assert 0 <= out[k] <= 1
    assert len(out["per_position_mismatch_full_length"]) == 6
    assert out["length_ratio"] == pytest.approx(out["mean_length"] / 6)


def test_core_cli_smoke(tmp_path, capsys):
    out = tmp_path / "r.json"
    core.main(["--target", "ACGT", "--cells", "2", "--procs", "1", "--terminators",
               "--set", "n_strands=40", "--set", "volume_fl=1e4", "--out", str(out)])
    data = json.loads(out.read_text())
    assert data["params"]["terminators"] is True and "exact_match" in data["results"]
    with pytest.raises(SystemExit):
        core.main(["--set", "not_a_param=1"])


# ============================ 3. statistics code ================================
def test_latin_hypercube_stratification():
    rng = np.random.default_rng(0)
    u = sens.latin_hypercube(50, 4, rng)
    assert u.shape == (50, 4) and u.min() >= 0 and u.max() < 1
    for j in range(4):
        assert sorted((u[:, j] * 50).astype(int)) == list(range(50))


def test_scale_unit_bounds_and_log_midpoint():
    u = np.array([[0.0] * 4, [0.5] * 4])
    x = sens.scale_unit(u)
    for j, nm in enumerate(sens.NAMES):
        lo, hi, sc = sens.PARAM_SPACE[nm]
        assert x[0, j] == pytest.approx(lo)
        assert x[1, j] == pytest.approx(np.sqrt(lo * hi) if sc == "log" else (lo + hi) / 2)


def test_rank_with_ties():
    assert sens.rank([10, 20, 20, 30]).tolist() == [1.0, 2.5, 2.5, 4.0]


def test_spearman_known_values():
    assert sens.spearman([1, 2, 3, 4, 5], [5, 6, 7, 8, 7]) == pytest.approx(0.8207826816681233)
    assert sens.spearman([1, 2, 3, 4], [1, 4, 9, 16]) == pytest.approx(1.0)
    assert sens.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    assert np.isnan(sens.spearman([1, 2, 3], [5, 5, 5]))


def test_spearman_matches_scipy_if_available():
    sp = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(1)
    x, y = rng.random(60), rng.random(60) + 0.5 * rng.random(60)
    assert sens.spearman(x, y) == pytest.approx(sp.spearmanr(x, y)[0])


def test_src_recovers_linear_model():
    rng = np.random.default_rng(2)
    X = rng.random((2000, 3))
    y = 3 * X[:, 0] - 1 * X[:, 1] + 0.0 * X[:, 2] + 0.01 * rng.standard_normal(2000)
    beta, r2 = sens.src(X, y)
    expect = np.array([3, -1, 0]) * X.std(axis=0) / y.std()
    assert np.allclose(beta, expect, atol=0.02) and r2 > 0.99


def test_first_order_index_known_functions():
    rng = np.random.default_rng(3)
    X = rng.random((6000, 3))
    y = X[:, 0] + 0.5 * X[:, 1]                       # true S1 = 0.8, 0.2, 0
    s = [sens.first_order_index(X[:, j], y) for j in range(3)]
    assert s[0] == pytest.approx(0.8, abs=0.04)
    assert s[1] == pytest.approx(0.2, abs=0.04)
    assert s[2] < 0.01


def test_s1_detects_non_monotone_effect_that_spearman_misses():
    rng = np.random.default_rng(4)
    x = rng.random(4000)
    y = (x - 0.5) ** 2
    assert abs(sens.spearman(x, y)) < 0.05
    assert sens.first_order_index(x, y) > 0.9


def test_global_indices_synthetic():
    rng = np.random.default_rng(5)
    U = sens.latin_hypercube(1500, 4, rng)
    Y = np.column_stack([U[:, 0] * 2 - U[:, 1], U[:, 2], np.zeros(1500) + 1e-9 * U[:, 3]])
    g = sens.global_indices(U, Y, n_boot=40, rng=rng)
    p = g["exact_match"]["params"]
    assert p["k_act"]["spearman"] > 0.8 and p["k_rev"]["spearman"] < -0.3
    assert g["edit_per_target_base"]["params"]["dark_leak"]["s1"] > 0.9
    lo, hi = p["k_act"]["spearman_ci95"]
    assert lo <= p["k_act"]["spearman"] <= hi


# ============================ 4. physics & end-to-end ===========================
def _metric(mode, metric, seed=0, **over):
    base = Params(n_strands=400, volume_fl=800.0, **sens.MODES[mode])
    out = analyse(simulate_cell(replace(base, **over), "ACGTAC", seed), "ACGTAC")
    return out[metric]


def test_bias_ratio_improves_fidelity_with_terminators():
    assert _metric("terminator", "exact_match", bias_ratio=200) > \
        _metric("terminator", "exact_match", bias_ratio=1) + 0.5


def test_more_light_gives_more_extension_without_terminators():
    lo = _metric("free", "length_ratio", k_act=0.2)
    hi = _metric("free", "length_ratio", k_act=10.0)
    assert hi > lo * 1.3


def test_dark_leak_increases_extension_without_terminators():
    assert _metric("free", "length_ratio", dark_leak=0.2) > _metric("free", "length_ratio", dark_leak=0.001) * 1.3


def test_faster_reversion_reduces_extension_without_terminators():
    assert _metric("free", "length_ratio", k_rev=5.0) < _metric("free", "length_ratio", k_rev=0.01)


def test_free_regime_baseline_is_not_capacity_saturated():
    p = Params(n_strands=300, **sens.MODES["free"])
    r = simulate_cell(p, "ACGTAC", 0)
    cap = int(6 * p.max_len_factor) + 1
    assert r["length"].mean() < 0.6 * cap and (r["length"] >= cap).mean() < 0.05


def test_evaluate_shapes_and_determinism():
    base = Params(n_strands=40, volume_fl=1e4, **sens.MODES["terminator"])
    pts = sens.scale_unit(sens.latin_hypercube(5, 4, np.random.default_rng(0)))
    m1, s1 = sens.evaluate(pts, sens.NAMES, base, "ACGT", 2, 3, 1)
    m2, _ = sens.evaluate(pts, sens.NAMES, base, "ACGT", 2, 3, 1)
    assert m1.shape == (5, 3) == s1.shape and np.array_equal(m1, m2)


def test_oat_structure_and_swing():
    base = Params(n_strands=40, volume_fl=1e4, **sens.MODES["free"])
    o = sens.oat_sweep(base, "ACGT", levels=3, reps=1, seed=1, procs=1)
    assert o["baseline"]["bias_ratio"] == base.bias_ratio
    for nm in sens.NAMES:
        e = o["params"][nm]
        lo, hi, _ = sens.PARAM_SPACE[nm]
        assert e["grid"][0] == pytest.approx(lo) and e["grid"][-1] == pytest.approx(hi)
        m = e["mean"]["length_ratio"]
        assert e["swing"]["length_ratio"] == pytest.approx(max(m) - min(m))


def test_run_sweep_terminator_regime_end_to_end(tmp_path):
    res = sens.run_sweep("terminator", n_strands=120, lhs=48, reps=1, levels=3, seed=3, n_boot=60)
    g = res["global_"]["exact_match"]["params"]
    top = max(g, key=lambda n: g[n]["s1"])
    assert top == "bias_ratio" and g["bias_ratio"]["spearman"] > 0.8
    md = sens.to_markdown(res)
    assert all(nm in md for nm in sens.NAMES)
    out = tmp_path / "s.json"
    sens.main(["--mode", "terminator", "--lhs", "16", "--reps", "1", "--levels", "3",
               "--n-strands", "40", "--procs", "1", "--json", str(out)])
    assert json.loads(out.read_text())[0]["mode"] == "terminator"


def test_run_sweep_free_regime_light_params_dominate():
    res = sens.run_sweep("free", n_strands=120, lhs=48, reps=1, levels=3, seed=4, n_boot=60)
    g = res["global_"]["length_ratio"]["params"]
    top = max(g, key=lambda n: g[n]["s1"])
    assert top in ("k_act", "dark_leak", "k_rev")
    assert g["k_act"]["spearman"] > 0 and g["dark_leak"]["spearman"] > 0 and g["k_rev"]["spearman"] < 0


def _light_share(res, resp="edit_per_target_base"):
    g = res["global_"][resp]["params"]
    return sum(g[n]["s1"] for n in ("k_act", "k_rev", "dark_leak")), g


def test_light_parameters_matter_only_when_light_dose_is_limiting():
    """Saturated terminator chemistry hides the light parameters (their combined S1 sits
    at the estimator noise floor, ~0-0.07 over seeds 5-7); slow terminator incorporation
    makes them clearly measurable (~0.27-0.51 over the same seeds)."""
    kw = dict(n_strands=150, lhs=64, reps=1, levels=3, seed=5, n_boot=40)
    sat, g_sat = _light_share(sens.run_sweep("terminator", **kw))
    slow, g_slow = _light_share(sens.run_sweep("terminator_slow", **kw))
    assert sat < 0.12 and slow > 0.20 and slow > sat + 0.15
    assert g_sat["bias_ratio"]["s1"] > 0.8
    assert g_slow["k_act"]["spearman"] < 0 < g_slow["k_rev"]["spearman"]


def test_run_sweep_rejects_unknown_mode():
    with pytest.raises(ValueError):
        sens.run_sweep("bogus")
```

## 7. Full test log

<details><summary>pytest -v (54 passed)</summary>

```text
test_tdt.py::test_params_validation_rejects[bad0] PASSED                
test_tdt.py::test_params_validation_rejects[bad1] PASSED                
test_tdt.py::test_params_validation_rejects[bad2] PASSED                
test_tdt.py::test_params_validation_rejects[bad3] PASSED                
test_tdt.py::test_params_validation_rejects[bad4] PASSED                
test_tdt.py::test_params_validation_rejects[bad5] PASSED                
test_tdt.py::test_params_validation_rejects[bad6] PASSED                
test_tdt.py::test_params_validation_rejects[bad7] PASSED                
test_tdt.py::test_params_validation_rejects[bad8] PASSED                
test_tdt.py::test_params_validation_rejects[bad9] PASSED                
test_tdt.py::test_params_validation_rejects[bad10] PASSED               
test_tdt.py::test_params_validation_rejects[bad11] PASSED               
test_tdt.py::test_params_validation_rejects[bad12] PASSED               
test_tdt.py::test_params_validation_rejects[bad13] PASSED               
test_tdt.py::test_params_validation_rejects[bad14] PASSED               
test_tdt.py::test_params_validation_rejects[bad15] PASSED               
test_tdt.py::test_params_validation_rejects[bad16] PASSED               
test_tdt.py::test_params_validation_rejects[bad17] PASSED               
test_tdt.py::test_params_defaults_valid_and_json_roundtrip PASSED       
test_tdt.py::test_bad_target_rejected PASSED                            
test_tdt.py::test_levenshtein_known_answers PASSED                      
test_tdt.py::test_truncated_product_counts_as_error PASSED              
test_tdt.py::test_max_run PASSED                                        
test_tdt.py::test_activity_step_matches_numeric_ode PASSED              
test_tdt.py::test_activity_step_zero_rates_is_identity PASSED           
test_tdt.py::test_poisson_limit_constant_activity PASSED                
test_tdt.py::test_mean_extension_matches_analytic_activity_integral PASSED
test_tdt.py::test_perfect_terminators_give_exact_product PASSED         
test_tdt.py::test_failed_deblock_creates_deletions_not_insertions PASSED
test_tdt.py::test_mass_conservation_under_depletion PASSED              
test_tdt.py::test_carryover_never_creates_molecules PASSED              
test_tdt.py::test_reproducible_and_seed_sensitive PASSED                
test_tdt.py::test_analyse_outputs_in_range PASSED                       
test_tdt.py::test_core_cli_smoke PASSED                                 
test_tdt.py::test_latin_hypercube_stratification PASSED                 
test_tdt.py::test_scale_unit_bounds_and_log_midpoint PASSED             
test_tdt.py::test_rank_with_ties PASSED                                 
test_tdt.py::test_spearman_known_values PASSED                          
test_tdt.py::test_spearman_matches_scipy_if_available PASSED            
test_tdt.py::test_src_recovers_linear_model PASSED                      
test_tdt.py::test_first_order_index_known_functions PASSED              
test_tdt.py::test_s1_detects_non_monotone_effect_that_spearman_misses PASSED
test_tdt.py::test_global_indices_synthetic PASSED                       
test_tdt.py::test_bias_ratio_improves_fidelity_with_terminators PASSED  
test_tdt.py::test_more_light_gives_more_extension_without_terminators PASSED
test_tdt.py::test_dark_leak_increases_extension_without_terminators PASSED
test_tdt.py::test_faster_reversion_reduces_extension_without_terminators PASSED
test_tdt.py::test_free_regime_baseline_is_not_capacity_saturated PASSED 
test_tdt.py::test_evaluate_shapes_and_determinism PASSED                
test_tdt.py::test_oat_structure_and_swing PASSED                        
test_tdt.py::test_run_sweep_terminator_regime_end_to_end PASSED         
test_tdt.py::test_run_sweep_free_regime_light_params_dominate PASSED    
test_tdt.py::test_light_parameters_matter_only_when_light_dose_is_limiting PASSED
test_tdt.py::test_run_sweep_rejects_unknown_mode PASSED                 
============================= 54 passed in 11.44s ==============================
```

</details>
