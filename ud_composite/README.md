# UD Composite Layup Sweep

Specimens **#291** and **#197** (가장 안전한 두 샘플) 의 동일한 cutout 형상을
유지하면서, woven fabric (`GeneralStiffnessSection`) 대신 **CU250 UD**
복합재 + `CompositeShellSection` 으로 적층각을 바꿔가며 Abaqus 동적해석을
수행한다.

## 디렉터리 구성

```
ud_composite/
  ud_layups.py            # 74개 적층각 카탈로그 (symmetric [...]s 전개됨)
  build_ud_specimen.py    # Abaqus CAE: 형상 + UD + 적층각 → INP
  post_ud_specimen.py     # Abaqus CAE: ODB → TSAIH/TSAIW/MSTRS max
  run_ud_all.py           # 드라이버: 적층각별 2-병렬 + 후처리
  sim_L<id>_S<num>/       # 자동 생성. INP/ODB/STA/MSG 등
  results_ud.csv          # 자동 생성. 적층각별 1행
```

## 해석 사양

| 항목 | 값 |
| --- | --- |
| 솔버 | Abaqus 2017 Explicit |
| Step1 (pinching) | 1.0 s |
| Step2 (folding)  | 3.0 s |
| 코어 / job | 6 |
| 동시 job | 2 (specimen 291, 197 병렬) |
| 형상 | `specimen/cutout_points_<num>.csv` (291, 197) |
| 재료 | CU250 UD (LAMINA + FailStress) |
| 단면 | `CompositeShellSection` (ply별 0.23 mm, 3 int pts, Simpson) |
| Material orientation | 원통 좌표계, ply 1축 = 글로벌 Z (축방향) |

## 재료 물성 (consistent units: t, mm, N, s, MPa)

| 파라미터 | 값 |
| --- | --- |
| ρ | 1.984e-9 t/mm³ |
| E₁ | 135 000 MPa |
| E₂ | 8 600 MPa |
| ν₁₂ | 0.30 |
| G₁₂ | 4 100 MPa |
| G₁₃ = G₂₃ | 2 000 MPa |
| Xₜ / X꜀ | 2 916 / 1 127 MPa |
| Yₜ / Y꜀ | 55 / 180 MPa |
| S | 113 MPa |
| ply 두께 | 0.23 mm |

## 파손기준

Abaqus 내장 `CFAILURE` 출력 (FailStress 정의 시 자동 계산) 의 다음
세 성분을 모든 frame · 모든 element · 모든 section point 에서 최댓값으로
수집한다.

| CFAILURE 성분 | 파손이론 |
| --- | --- |
| `TSAIH` | Tsai–Hill |
| `TSAIW` | Tsai–Wu (F₁₂ = -0.5√(F₁₁·F₂₂)) |
| `MSTRS` | Maximum Stress |

값 ≥ 1.0 이면 해당 ply가 파손했다는 의미. 각 layup ID 당 한 행을 갱신한다.

## 실행 예

```powershell
# 전체 74 적층각 (12-core PC에서 6×2 병렬)
abaqus python ud_composite\run_ud_all.py

# 4-ply 그룹만 (id 36..54)
abaqus python ud_composite\run_ud_all.py --layup-start 36 --layup-end 54

# 특정 ID들만
abaqus python ud_composite\run_ud_all.py --layup-ids 11,36,55,69

# 291만 (197 패스)
abaqus python ud_composite\run_ud_all.py --specimens 291
```

각 적층각의 진행 순서:

1. **Phase 1** — specimen 291 INP build → 197 INP build (순차)
2. **Phase 2** — solver `ud_L<id>_S291` + `ud_L<id>_S197` 병렬 실행
3. **Phase 3** — 두 ODB 순차 후처리, `results_ud.csv` 동일 행 업데이트

## results_ud.csv 스키마

| 컬럼 | 설명 |
| --- | --- |
| LAYUP_ID | 1..74 |
| GROUP | `2-ply` / `3-ply` / `4-ply` / `6-ply` / `8-ply` |
| LAYUP_LABEL | `[45/-45]s` 등 원본 표기 |
| PLY_COUNT | 정수 |
| LAYUP_ANGLES | `45/-45/-45/45` 형식 (실제 사용 각도) |
| TSAI_HILL_291, TSAI_HILL_197 | 두 specimen의 Tsai-Hill max |
| TSAI_WU_291,   TSAI_WU_197   | Tsai-Wu max |
| MAX_STRESS_291, MAX_STRESS_197 | Max-stress 기준 max |

## 참고

- `build_ud_specimen.py` 의 형상/BC/하중/메시 로직은
  `../tube_hinge_pinching.py` 와 동일. 차이점은 [6] PROPERTY 모듈
  (woven → UD CompositeShellSection) 과 field output에 `CFAILURE`
  추가뿐.
- `post_ud_specimen.py` 는 기존 FI 계산 대신 Abaqus 내장 CFAILURE
  값을 그대로 사용. component label로 TSAIH/TSAIW/MSTRS 인덱스를
  찾아 모든 frame에서 max를 누적한다.
- 빌드/포스트 스크립트 모두 Abaqus 2017 (Python 2.7) 호환을
  유지했다 (`execfile`, `print_function`, `wb` fallback 등).
