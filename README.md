# Tape Hinge — LHS Sampling

복합재 테이프 스프링 힌지의 Bezier C2 cutout 형상 40개에 대한 아바쿠스 해석 자동화 파이프라인.
MBO 루프 진입 전, surrogate 학습용 초기 데이터를 LHS로 추출하는 단계.

## 사양

| 항목 | 값 |
|---|---|
| 힌지 길이 `L_total` | 220 mm |
| 원통 직경 `D` | 38 mm |
| 면밀도 | 318 g/m² (`3.18e-10` kg/mm²) |
| 샘플 수 | 40 (Bezier C2 6CP, LHS) |
| 설계변수 | `Eccentricity`, `y1`, `x1`, `Dx`, `y3`, `yP2`, `fP3` |
| 파생변수 | `x2_derived`, `x3_derived`, `xP3_derived`, `yP3_derived` |
| 해석 | Pinching (1s) + Folding (3s), Explicit Dynamics |

## 빠른 시작

```bash
git clone <repo-url> Tape_hinge_MBO_LHS_Sampling
cd Tape_hinge_MBO_LHS_Sampling
python run_all.py --start <S> --end <E> --batch <B> --cpus <C>
```

스크립트는 자기 위치(`__file__`)에서 모든 경로를 도출하므로 어디서 clone해도 추가 설정 불필요.

## 컴퓨터별 실행 명령

40개 샘플을 3대로 분산:

| PC | 코어/job | 동시 job | 샘플 | 명령 |
|---|---|---|---|---|
| **PC1** | 6 | 3 | 0..18 (19개) | `python run_all.py --start 0  --end 18 --batch 3 --cpus 6` |
| **PC2** | 6 | 2 | 19..30 (12개) | `python run_all.py --start 19 --end 30 --batch 2 --cpus 6` |
| **PC3** | 6 | 2 | 31..39 (9개) | `python run_all.py --start 31 --end 39 --batch 2 --cpus 6` |

각 PC는 자신이 담당한 specimen만 처리해 로컬 `results.csv`에 기록.
3대 종료 후 사용자가 수동으로 `results.csv` 3개를 병합.

## 파일 구성

| 파일 | 역할 |
|---|---|
| `run_all.py` | 메인 배치 자동화 (Python 3, argparse) |
| `tube_hinge_pinching.py` | Abaqus CAE — Part/Property/Assembly/Step/Mesh/Job 통합 |
| `post_process_csv_ver2.py` | Abaqus CAE — ODB 후처리 → `results.csv` 누적 |
| `get_FI_graph_pinching_ver2.py` | Abaqus CAE — Failure Index 시각화 (수동) |
| `specimen/cutout_points_<N>.csv` | 40개 cutout 좌표 + DV |
| `bezier_c2_cutout_sampling.ipynb` | 샘플 생성 노트북 (참고) |
| `01_part.py` ~ `05_mesh.py` | 모듈형 레거시 (디버깅 참고용) |

## 파이프라인 단계 (batch별)

```
Phase 1: tube_hinge_pinching.py 로 INP 생성 (sequential, abaqus cae noGUI)
Phase 2: abaqus job=... 솔버 병렬 실행 (--batch개 동시)
Phase 3: post_process_csv_ver2.py 로 ODB → results.csv (sequential)
```

## 산출물

- `sim_<N>/cutout_points_<N>.odb` — 해석 결과 (~1.5 GB each, gitignored)
- `results.csv` — 누적 결과 (gitignored)
  - 컬럼: `SPECIMEN_NUM, MAX_FI_ALL, MAX_FI_PINCHING, MAX_FI_FOLDING, MAX_SE, Eccentricity, y1, x1, Dx, y3, yP2, fP3, x2_derived, x3_derived, xP3_derived, yP3_derived`

## 사전 조건

- Windows + Abaqus CAE 6.14+ (또는 호환 버전)
- `abaqus` 명령이 PATH에 등록되어 있어야 함
- Python 3 (drivers)

## 다음 단계

3대의 `results.csv`를 합쳐 40행 완성 → MBO 루프 레포지토리로 이동하여 surrogate 학습 + Bayesian Optimization 진행.
