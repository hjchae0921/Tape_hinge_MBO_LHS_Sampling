# new_simulation — Folded-FI MBO (qLogEHVI)

A **self-contained** multi-objective Bayesian optimization campaign for the
tape-hinge Bezier cutout, identical to the parent `mbo/` pipeline **except the
failure objective is the FULLY-FOLDED failure index** instead of the
max-over-time FI.

## Objectives (both maximized)

| | |
|---|---|
| `y0 = 1 − MAX_FI_FOLDED` | `MAX_FI_FOLDED` = max-over-elements FI at the **last folding frame** (fully folded) |
| `y1 = MAX_SE` | max whole-model strain energy (ALLSE) |

Reference point `(0.0, 400.0)`, Matern12 ARD GP, fixed noise 1e-6, seed 42,
qLogEHVI, total budget 400 = 40 initial LHS + 360 BO.

`results.csv` also stores `MAX_FI_ALL / PINCHING / FOLDING` (max-over-time) for
reference, but the **optimizer only uses `MAX_FI_FOLDED`**.

## Layout (everything relative to this folder)

```
new_simulation/
  run_campaign.py          master: initial LHS 0..39 -> BO 40..399
  run_mbo.py               qLogEHVI BO loop (folded objective)
  gp_model.py              GP + qLogEHVI  (host python / botorch)
  bezier_cutout.py         cutout geometry
  abaqus_runner.py         3-phase Abaqus wrappers (Py2/Py3 io.open)
  plotter.py               Pareto scatter
  tube_hinge_pinching.py   Abaqus build  (Abaqus 2017 Py2 + 2024 Py3)
  post_process_folded.py   Abaqus post   (folded FI + SE)
  seed_existing.py         (this machine only) seed 0..11 from retained ODBs
  specimen/cutout_points_0..39.csv   shared initial 40 LHS geometries
  results.csv, sim_<N>/, plot/        generated (gitignored)
```

## Environment

**Host python** (runs the BO driver, NOT Abaqus):
```
conda create -n mbo python=3.11 -y
conda activate mbo
pip install -r requirements.txt          # torch, botorch, gpytorch, numpy, scipy, pandas
```
`abaqus` must be on PATH (`abaqus cae`, `abaqus job=...`).

## Run

**q=1 machine (Abaqus 2017):**
```
python new_simulation/run_campaign.py --q 1 --cpus 8
```
**q=3 machine (Abaqus 2024, 24 cores):**
```
python new_simulation/run_campaign.py --q 3 --cpus 8
```
Resumable — just re-run; finished designs in `results.csv` are skipped.
A fresh machine evaluates all 40 initial LHS first, then BO to 400.

## Abaqus 2017 (Python 2.7) compatibility notes

The Abaqus-side scripts (`tube_hinge_pinching.py`, `post_process_folded.py`)
and the generated wrappers are written for both Py2.7 (2017) and Py3 (2024):
`from __future__ import print_function`, `io.open(..., encoding=...)`, and a
`csv` `newline=''`→`'wb'` fallback.

⚠️ **Test one build first on the 2017 box** before the full run:
```
python new_simulation/run_campaign.py --q 1 --cpus 8   # watch the first build
```
If the build errors, the most likely 2017-vs-2024 API difference is the
`ExplicitDynamicsStep(..., massScaling=(( ... ),))` tuple arity — adjust the two
`massScaling=` tuples in `tube_hinge_pinching.py` to the 2017 signature if needed.
