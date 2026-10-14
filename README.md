# CPD-RNG: Compound Poisson Distribution Random Number Generator

Reference implementation for the manuscript:

> **Compound Poisson Modeling of System Timing Jitter for Statistically
> Certified Random Number Generation**
> Randy Kuang, Quantropi Inc.
> *Academia Quantum*

---

## Overview

This repository contains the Java reference implementation of the CPD
random number generator described in the manuscript. The construction
emulates the photon-number statistics of coherent states through a
Poisson-modulated random permutation process driven by physical timing
jitter. Each cycle draws a Poisson-distributed permutation count, executes
that many Fisher–Yates permutations, measures the cumulative elapsed
time, and advances a 48-bit LCG state by the modular residue of the
elapsed time. The emitted symbol is the permutation count modulo R.

The code reproduces every experiment reported in Section 3 of the
manuscript, including the Poisson sampler verification (Tables 2–3,
Figure 2), the CPD characterization under different system loads
(Figures 3–4), the same-seed time-series comparison (Tables 4–7),
the SUP validation (Table 8), and the uniformity and NIST SP 800-90B
analyses (Tables 10–14).

---

## Files

- `CPD_RNG.java` — the reference implementation
- `README.md` — this file
- `LICENSE` — license under which the code is distributed

---

## Requirements

- OpenJDK 17 or later
- No external dependencies

---

## Default configuration

The default parameters of the shipped code correspond to the
**high-count configuration** of the manuscript:

| Parameter     | Value            | Role                                  |
|---------------|------------------|---------------------------------------|
| `MU`          | `100.0`          | Poisson mean of the permutation count |
| `R`           | `16`             | output modulus (4-bit symbol)         |
| symbol bits   | `4`              | bits per symbol, passed to sampler    |
| `INIT_SEED`   | `347563409L`     | initial LCG seed                      |
| `INIT_ARRAY`  | `{3, 1, 2, 4}`   | initial array state                   |
| `WARMUP`      | `8`              | discarded cycles before recording     |

This configuration is the one used for the results reported in
Tables 10, 11, 13, and 14, and in Figures 2, 3, and 4 of the
manuscript.

---

## Switching to the low-count configuration (μ = 7, R = 4)

The manuscript also reports results for the low-count configuration
(`μ = 7`, `R = 4`), shown in Table 10 and in Figure 1. Running this
configuration requires changing **three parameters** in the source
code, because the code contains multiple references to `μ` and to
the symbol width that must remain consistent.

Open `CPD_RNG.java` and update the following three locations:

**1. The `MU` constant** (top of the class):

```java
public static final double MU = 7.0;   // was 100.0