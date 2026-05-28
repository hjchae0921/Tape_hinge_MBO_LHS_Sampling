# MBO — Multiobjective Bayesian Optimization for Tape Hinge

Initial 40 LHS samples are already in `../results.csv`. This driver runs the
remaining **200 BO iterations** (`SPECIMEN_NUM` 40 → 239) with qLogEHVI.

## Specification (from `../MBO.md`)

| 항목 | 값 |
|---|---|
| 총 예산 | 240 (초기 LHS 40 + BO 200) |
| 목적함수 (maximize) | `y0 = 1 − MAX_FI_ALL`, `y1 = MAX_SE` |
| Reference point | `(0.0, 400.0)` |
| Kernel | Matern12 (ν = 0.5), ARD |
| Observation noise | fixed `σ² = 1e-6` (deterministic; jitter for Cholesky) |
| Hyperparameters | 매 iteration **MAP** refit (BoTorch default GammaPrior — lengthscale Γ(3, 6), outputscale Γ(2, 0.15)) |
| Acquisition | **qLogEHVI** (qEHVI의 log-domain numerical fix variant — BoTorch 권장, arXiv:2310.20708. 수학적으로 동일), **q = 3** (joint batch) |
| Random seed | 42 |
| Solver | Abaqus 2024 Explicit, **4 cpus/job × 3 jobs 병렬** (= 12 cores 활용) |
| Build / Post | sequential (CAE license · `results.csv` 동시쓰기 회피) |

## 설치

MBO 드라이버는 **호스트 Python** 에서 돌고, Abaqus 자식 프로세스만 호출합니다.
Abaqus 내부 인터프리터와 분리된 환경을 만드세요.

**권장: conda**

```powershell
cd C:\Users\user\Desktop\Tape_hinge_MBO_LHS_Sampling
conda create -n mbo python=3.11 -y
conda activate mbo
pip install -r mbo\requirements.txt
```

GPU(CUDA)를 쓰실 거면 `pip install -r` 이전에 PyTorch를 먼저 설치:

```powershell
conda install pytorch pytorch-cuda=12.1 -c pytorch -c nvidia -y
pip install -r mbo\requirements.txt
```

(BO overhead는 솔버 시간 대비 작아 CPU만으로도 충분합니다.)

**대안: venv**

```powershell
python -m venv .venv-mbo
.\.venv-mbo\Scripts\activate
pip install -r mbo\requirements.txt
```

`abaqus` 명령이 PATH에 있어야 합니다 (`abaqus cae`, `abaqus job=...`).
이후 매번 새 PowerShell에서 작업할 때 환경 활성화부터:

```powershell
conda activate mbo            # 또는: .\.venv-mbo\Scripts\activate
python mbo\run_mbo.py
```

## 실행

```powershell
python mbo\run_mbo.py                    # 200개 specimen 자동 진행 (q=3, cpus=4)
python mbo\run_mbo.py --max-iters 6      # smoke test (specimen 6개만 → 2 batch)
python mbo\run_mbo.py --q 3 --cpus 4     # 명시적 설정
python mbo\run_mbo.py --q 2 --cpus 6     # 다른 PC 사양에 맞춰 변경
```

마지막 batch는 남은 budget에 맞춰 q를 자동으로 줄입니다 (200 / 3 → 66 batch + 마지막 q=2).

## 출력

| 파일 | 내용 |
|---|---|
| `../results.csv` | 후처리 스크립트가 매 iteration 후 갱신 (16열) |
| `../mbo_log.csv` | iteration별 status, HV, qLogEHVI, DV, GP hyperparameters |
| `../plot/pareto_latest.png` | 가장 최근 batch 종료 직후 scatter + Pareto (자동 덮어쓰기) |
| `../plot/pareto_iter_<NNN>.png` | batch 단위 history snapshot |
| `../sim_<N>/` | Abaqus 작업 파일 (gitignored) |
| `../specimen/cutout_points_<N>.csv` | BO가 제안한 cutout 좌표 + DV |

`plot/` 은 `.gitignore`로 제외되어 있어 커밋되지 않습니다. scatter는 x = `1 − MAX_FI_ALL`,
y = `MAX_SE` 의 maximize 평면으로 그리고, 초기 LHS는 회색, BO 점은 파란색, Pareto front는
빨간색 빈 동그라미 + 단계선, reference (0, 400) 은 검정 별표로 표시됩니다.

## 동작 흐름 (한 batch iteration, q=3 기준)

```
1. results.csv 다시 읽기 → (X, Y) tensor + next_idx = max(SPECIMEN_NUM)+1
2. GP fit (2개 SingleTaskGP, fixed noise σ²=1e-6, Matern12, **MAP** with default GammaPriors)
3. qLogEHVI optimize_acqf → 3 × 7-DV 후보 joint 제안
4. bezier_cutout.write_cutout_csv ×3   → specimen/cutout_points_{N,N+1,N+2}.csv
5. [phase 1, 순차]   abaqus cae noGUI=<tube_hinge_pinching wrapper>  ×3
6. [phase 2, 병렬]   abaqus job=cutout_points_<N> cpus=4 interactive  ×3 동시
7. [phase 3, 순차]   abaqus cae noGUI=<post_process_csv_ver3 wrapper> ×3
8. 실패한 specimen 만 개별로 retry (sequential, --max-retries 회)
9. mbo_log.csv 에 3행 기록
```

## 실패 처리

- batch 단계에서 어떤 specimen이 phase 1/2/3 중 rc != 0 또는 ODB 누락 → 그 specimen만
  표시 후 phase 진행을 멈춤 (다른 specimen은 계속).
- batch 종료 후 실패 specimen은 sim_<N>/ 폴더를 정리하고 **개별 sequential retry 1회**.
- retry까지 실패 → `results.csv`에 placeholder 행 (`SPECIMEN_NUM`만 채움) 추가하여
  다음 batch가 같은 index에 갇히지 않도록 함. `mbo_log.csv`에 `status=failed`.

## 재개 (resume)

`Ctrl+C`로 중단하거나 PC가 꺼져도, 다시 `python mbo\run_mbo.py` 만 호출하면
`results.csv` 마지막 `SPECIMEN_NUM` + 1 부터 자동 재개합니다.

## 주의

- `sim_<N>` 또는 `specimen/cutout_points_<N>.csv` 가 이미 존재하면 driver는
  덮어쓰지 않고 abort합니다. (의도치 않은 손실 방지) 재실행 전에 수동으로 정리하세요.
- `results.csv` 의 trailing empty 컬럼 (현재 3개) 는 첫 번째 BO iteration의
  post-process 단계에서 자동으로 16-열 표준 헤더로 정리됩니다.
