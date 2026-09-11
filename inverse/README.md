# High-NA scalar ODT inverse experiment

Three reconstructions of the same independently generated Maxwell
co-polarized detector fields: direct Fourier gridding, nonnegative
Gerchberg–Papoulis (GP), and nonnegative TV with Lim residual reinjection.
The selected system is the completed forward study's Δn × 0.01 case.
See [results/report.md](results/report.md) for measured errors and limitations;
[high_na_plan.md](high_na_plan.md) records the approved scope.

## 현재까지의 구성과 도달점 — 2026-09-13

**Δn × 0.01 조건의 11개 run·17개 복원과 별도의 sampling 검사를 완료했다.**
이는 예정한 계산이 끝났다는 뜻이며, 모든 알고리즘의 수렴이나 정확한
굴절률 복원을 달성했다는 뜻은 아니다. 아래는 현재 상태의 요약이고,
개별 수치·단면·원인 분석은 [최종 보고서](results/report.md)에 있다.

### 풀고 있는 문제와 구현 구성

대상은 a = 3 μm(z), b = 5 μm(x, y)의 oblate spheroid이며,
λ₀ = 0.532 μm, nₘ = 1.335381534, nₚ = 1.33567771866,
Δn = 0.00029618466이다. 조명 NA 0.5·검출 NA 1.0에서 81개 입사 방향을
사용하고, 빛의 진행 방향은 kᵢ,z > 0, 검출면은 고정된 z = +5 μm다.
이 Δn는 원래 contrast를 인위적으로 0.01배 낮춘 검증용 조건이다.

**a, b, nₚ 세 변수만 추정하는 fitting이 아니라, 128 × 128 × 80개 복셀의
3D 굴절률 분포를 복원한다.** 복셀 간격은 0.1 μm이며 검출 격자는 별도의
512 × 512·0.1 μm다. 실제 spheroid support, 균질성, 굴절률 상한은
inverse에 주지 않았고, 반복 방법에는 비음수 contrast 조건 χ ≥ 0을 준다.
잡음은 추가하지 않았으며, 모든 주 실험은 Rytov 전처리 데이터를 사용한다.

| 구성 파일 | 현재 역할 |
|---|---|
| [odt_exact_data.m](odt_exact_data.m) | 완료된 exact Maxwell reference와 기존 편광·co-projection을 재사용하여 독립 관측값 구성 |
| [odt_prepare.m](odt_prepare.m) | 실제 pupil을 통과한 total field의 Rytov logarithm과 Ewald 샘플 생성 |
| [odt_operator.m](odt_operator.m) | Padded Fourier transform·trilinear interpolation 및 정확히 대응하는 adjoint |
| [odt_reconstruct.m](odt_reconstruct.m) | Direct Fourier, 비음수 GP, Lim TV의 세 복원 방법 |
| [odt_predict.m](odt_predict.m) | 복원값을 기존 Born/Rytov 물리식으로 다시 forward 계산하여 산란장 오차 평가 |
| [run_inverse_experiment.m](run_inverse_experiment.m) | 설정 검사, checkpoint 재사용, 복원 및 오차·이력 저장 |
| [check_inverse_sampling.m](check_inverse_sampling.m) | Fourier padding, 복셀 간격, 검출 pitch/window 검사 |
| [plot_inverse_results.py](plot_inverse_results.py) | 저장된 RI 오차의 독립 재계산과 단면·이력·missing cone 그림 생성 |

Lim TV는 논문의 식 (14), (15)와 Appendix A의 splitting을 따르되,
현재 비균일 Ewald 연산자에 맞는 선형계를 PCG로 풀고 수렴 조건을 검사한다.
정규화 강도는 관측 데이터에서 얻는 크기에 비례하도록 설정했으며,
정답 형상이나 정답 Δn로 최적화하지 않았다.

### 수행한 실험: 세 알고리즘, 11개 설정, 17개 복원

| Run 이름 | 방법 / 복원 수 | 목적 |
|---|---|---|
| `exact` | Direct·GP·TV / 3 | 독립 exact 관측값, padding 4의 기본 비교 |
| `linear_control` | Direct·GP·TV / 3 | 같은 이산 연산자로 y = Aχ_true를 합성하여 물리적 모델 불일치를 제거한 대조군 |
| `padding8` | Direct·GP·TV / 3 | 기본 반복 한도를 유지하며 Fourier interpolation 정밀도 증가 |
| `gp_control_tight` | GP / 1 | Matched 대조군에서 최대 300회, 더 엄격한 CGLS 계산 |
| `gp_tight` | GP / 1 | Exact 데이터에서 동일한 엄격한 GP 설정 |
| `tv_long` | TV / 1 | 기본 TV 강도, outer당 최대 1000회로 inner 반복 증가 |
| `tv_half` | TV / 1 | 긴 inner 계산에서 TV 강도를 기본의 ½로 변경 |
| `tv_double` | TV / 1 | 긴 inner 계산에서 TV 강도를 기본의 2배로 변경 |
| `tv_penalty` | TV / 1 | 같은 정규화 목적함수에서 splitting penalty ρ = ρ₊ = 0.01 검사 |
| `tv_refined` | TV / 1 | Padding 8에서 data·TV·splitting penalty의 상대적 크기를 맞춘 검사 |
| `tv_strict` | TV / 1 | 같은 refined 문제에서 수렴 허용오차를 10⁻⁵로 강화 |
| **합계** | **17개 복원** | |

이 17개는 물체·Δn·조명/검출 NA를 바꾸지 않았다. Matched control의
정답 사용은 관측값 합성에 한정되며, inverse의 제약으로 넣은 것이 아니다.
추가 sampling 검사는 padding 2/3/4/6/8, 복셀 간격 0.1/0.08 μm,
대표 7개 조명의 detector pitch/window 변경이며, 별도의 복원 17개에
합산하지 않는다. 초기 구현과 중단된 시도는
[provisional 기록](results/provisional/README.md)으로 분리했다.

### 현재 결과와 결론

오차는 ε_n = ‖n_rec − n_true‖₂ / ‖n_true − nₘ‖₂로 정의하며,
큰 배경 굴절률이 아니라 실제 contrast로 정규화한 전체 부피 오차다.
아래는 같은 padding 8 연산자에서의 비교이며, GP는 기본 100회 한도,
TV는 penalty를 맞추고 tolerance 10⁻⁵를 적용한 `tv_strict` 결과다.

| 방법 | RI 오차 | 물리적 산란장 오차 | 수치 상태 |
|---|---:|---:|---|
| Direct Fourier | 46.9056% | 17.6344% | 직접 변환 |
| 비음수 GP | 36.3767% | 0.2263% | 100회 한도, 미수렴 |
| Lim TV, strict | 28.9450% | 0.9293% | 856 inner / 1 outer, 수렴 |

- TV가 이 비교에서 가장 정확했지만, 평균 Δn가 15.2254% 낮고 중심선의
  half-contrast 영역이 축방향 계산 경계에 닿았다. 보고된 8 μm span은
  정확히 복원된 물체 두께가 아니라 경계에 제한된 값이다.
- TV의 tolerance를 10⁻⁴에서 10⁻⁵로 강화해도 RI 오차는
  29.3580% → 28.9450%로만 변했다. 현재 약 29% 오차를 단순히 inner
  반복 부족만으로 설명할 수는 없다.
- Padding 4의 `tv_penalty`는 19.5501%였지만 padding refinement에서
  유지되지 않았다. 이를 정답 오차가 가장 작다는 이유로 최종 검증
  성능으로 선택하지 않았다.
- `gp_tight`는 300회에서 미수렴, RI 오차 85.0047%였으며 내부 CGLS
  correction 189회가 요구 정확도를 충족하지 못했다. 같은 엄격한 설정의
  matched control은 수렴했지만 RI 오차가 33.0805% 남았다.
- 최종 17개 중 Direct 3개를 제외한 반복 복원 14개는 수치 기준 충족
  7개·미충족 7개다. `INVERSE_EXPERIMENT_PASS`는 결과 저장과 실행 완료
  표시이며, 모든 복원의 수렴을 뜻하지 않는다.

**관측장을 잘 맞추는 것과 굴절률을 정확히 복원하는 것은 다르다.**
Missing cone에 따른 정보 부족, 모델·데이터 불일치, 이산화·보간,
GP의 불완전한 projection, TV 편향과 경계 조건을 함께 고려해야 한다.
현재 오차 전체를 scalar approximation이나 missing cone 하나로
귀속시키는 결론은 검증되지 않았다.

### Padding 2의 6.33%, padding 8의 0.63%가 뜻하는 것

Padding은 χ 배열 바깥에 0을 붙여 Fourier 계산점을 촘촘하게 만드는
계산 방법이며, 원래 물체·복셀 간격·관측 NA는 바꾸지 않는다.
한 축의 Fourier 격자 간격은 Δq = 2π/(pNΔx)이므로, padding 배수 p를
2에서 8로 늘리면 계산점 간격은 ¼이 되어 관측 좌표 사이의 보간이 개선된다.

정답 spheroid의 해석적 Fourier 값과 이산 계산값을 비교한 상대 오차가
6.3303% → 0.6314%로 줄었다는 뜻이며, **굴절률 복원 오차가 아니다.**
이 수치에는 순수한 보간 오차뿐 아니라 경계를 복셀로 표현한 오차도
포함된다. Padding은 새로운 측정 정보를 추가하거나 missing cone을
채우지 않으며, 이 검사가 RI 복원의 grid convergence를 증명하지도 않는다.
세부 값은 [operator_convergence.csv](results/operator_convergence.csv)에 있다.

### 검증 기록, 생성물, 아직 확인하지 못한 점

실험 종료 시점의 [최종 inverse 검증 로그](results/verification_final.log)는
회귀시험 PASS와 MATLAB Code Analyzer 12개 파일·0건을 기록한다.
[Artifact 검사](results/artifact_check.log)는 17개 복원의 MAT/CSV 설정,
이력, 수렴 flag와 독립 RI 오차 재계산 및 11개 그림을 확인한 기록이다.
별도의 기존 forward suite는 `test_born_rytov.m:178`의 tolerance 유효성
검사 불일치로 실패했으며, 그 함수/시험은 수정하지 않았다.

- `results/exact_detector_512.mat`: 주 실험이 공유하는 독립 관측 cache.
- `results/<run>_geometry.mat`: 격자·설정·평가용 truth·주파수 coverage.
- `results/<run>_<method>.mat`: 복원값, 설정, 오차, 수렴 진단.
- `results/<run>_metrics.csv`, `*_history.csv`, `<run>.log`: 수치 결과와 이력.
- `results/*.png`: 단면, 중심선, 반복 이력, missing cone 및 TV 민감도 그림.
- [results/report.md](results/report.md): 전체 결과와 해석;
  [verification.md](verification.md): 현재·과거 검증의 구분.

Periodic TV 경계의 기여, 더 큰 복원 상자, 정당한 background/support
조건의 효과, 실제 잡음에서의 안정성은 아직 원인별로 분리 검증하지 않았다.
TV 강도와 1% Ewald discrepancy는 선언한 설정이지 실제 장비 잡음에서
검증된 최적값이 아니다. 선택한 편광·co-projection의 정확도가 임의 편광과
실제 objective의 vector polarization transfer까지 보장하지도 않는다.
현재 실험은 완료했지만, 정확한 RI 복원과 이러한 원인의 정량적 분리는
미해결 상태다. 아래의 기존 기술 설명과 재현 명령을 그대로 유지한다.

## Conditions and reproduction

| Quantity | Value |
|---|---|
| Vacuum wavelength | 0.532 μm |
| Water / cell RI | 1.335381534 / 1.33567771866 |
| RI difference | 0.00029618466 |
| Oblate semiaxes | a = 3 μm along z, b = 5 μm along x/y |
| Illumination / detector NA | 0.5 / 1.0 |
| Illumination | 81 directions, positive kᵢ,z |
| Detector | fixed z = 5 μm, 512 × 512, pitch 0.1 μm |
| Reconstruction | 128 × 128 × 80, pitch 0.1 μm |
| Noise / object prior | none added / real nonnegative contrast only |

From the repository root:

```matlab
restoredefaultpath
addpath('inverse/tests')
run_all
addpath('inverse')
settings = struct('padding',4,'gp_iter',100,'tv_inner',200,'cg_max_iter',100);
run_inverse_experiment('exact',settings);
settings.source = 'matched';
run_inverse_experiment('linear_control',settings);
check_inverse_sampling;
```

The runner saves MAT checkpoints, CSV metrics/histories and a diary under
`results/`. Reusing a label requires identical settings; completed methods
are reused, and incomplete methods restart. The matched control uses Aχ_true
as its data and is explicitly **not** an independent physics validation.
Truth enters only this labeled control and evaluation, never a reconstruction
constraint or the primary TV parameter selection.
`shasum -a 256 -c inverse/results/source_data.sha256` checks the recorded
inverse source, reused forward helpers, selected reference and measurement cache
from the repository root.

To reproduce the complete 17-volume study, including its explicitly labeled
controls and unsuccessful numerical settings:

```matlab
base = struct('padding',4,'gp_iter',100,'tv_inner',200,'cg_max_iter',100,'rho',.1);
run_inverse_experiment('exact',base);
s = base; s.source = 'matched'; run_inverse_experiment('linear_control',s);
s = base; s.padding = 8; run_inverse_experiment('padding8',s);
g = base; g.methods = {'gp'}; g.gp_iter = 300;
g.cg_tol = 1e-5; g.cg_max_iter = 200;
g.source = 'matched'; run_inverse_experiment('gp_control_tight',g);
g.source = 'exact'; run_inverse_experiment('gp_tight',g);
t = base; t.methods = {'tv'}; t.tv_inner = 1000;
labels = {'tv_long','tv_half','tv_double'}; multipliers = [4,2,8];
for j = 1:3
    t.alpha_relative = multipliers(j);
    run_inverse_experiment(labels{j},t);
end
t.alpha_relative = 4; t.rho = .01;
run_inverse_experiment('tv_penalty',t);
a = load('inverse/results/exact_geometry.mat','physical_normalization');
b = load('inverse/results/padding8_geometry.mat','physical_normalization');
t.padding = 8; t.rho = .01*b.physical_normalization/a.physical_normalization;
run_inverse_experiment('tv_refined',t);
t.tol = 1e-5; t.cg_tol = 1e-5;
run_inverse_experiment('tv_strict',t);
```

The no-options runner is a coarse padding-2 baseline, not the refined setting.
Changing a run label's parameters is rejected before its artifacts are touched.

Plot and independently recompute saved RI metrics with the already available
NumPy/SciPy/Matplotlib environment:

```sh
/tmp/mie-oblate-plot-env/bin/python inverse/plot_inverse_results.py
```

The commands above restore the path for that process only; no preferences or
global path are saved. A no-desktop invocation is:

```sh
/Applications/MATLAB_R2024a.app/bin/matlab -nojvm -nodesktop -nosplash -r "restoredefaultpath; addpath('inverse/tests'); try, run_all; catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0)"
```

## Data and forward model

Array order is (x,y,z), with the origin at the grid center. Light arrives from
the −z side and has kᵢ,z > 0; tilted illumination does not rotate the detector.
NA = nₘ sin θ, so the selected half-angles are 21.989° and 48.491°.
The detector pupil is |kₛ,⊥| ≤ k₀ NA_det; object frequencies are q = kₛ − kᵢ.

Definitions and transform conventions:

- χ = n²/nₘ² − 1, f = kₘ²χ, n = nₘ√(1 + χ).
- Time dependence exp(−iωt); outgoing wave exp(+ikₘr).
- χ̂(q) = ∫ χ(r) exp(−iq·r) d³r.
- L_B = (U − u₀)/a; L_R = (u₀/a) log(U/u₀), where a is the reference amplitude.
- g(q) = L̂(kₛ,⊥)(−2ikₛ,z) exp(−ikₛ,z z_det)/kₘ².

`odt_exact_data` reuses coefficients and polarization conventions from
`../spheroid-analytic-forward/2/case_3_reference.mat` via
`../forward/exp/oblate_na_sweep.m`. It checks the completed reference and
parameter/angle identity, then constructs physical pupil-filtered detector
fields independently of the voxel inverse.
`odt_prepare` takes their total-field Rytov logarithm after the pupil, retains
the nonuniform Ewald samples, and checks reference amplitude and phase branches.
Automatic phase unwrapping compares x/y path orders and requires a near-zero
background border; it is a consistency check, not a general vortex/sampling proof.

`odt_predict` uses the same physical equations as `born_rytov_voxel`, with
independent volume and detector sizes: compute the pre-pupil Born field B,
ψ = B/u₀, then apply the pupil to u₀ expm1(ψ).
Filtering and logarithms/exponentiation do not commute, so physical post-pupil
Rytov data are not an exactly linear function of χ.
The pre-pupil hemisphere requires detector pitch < λ₀/(2nₘ) = 0.199194 μm.

The fixed transverse polarization and co-projection match the validated
forward study; the scalar inverse itself has no polarization argument.
Transversality E₀·kᵢ = 0 still applies. Accuracy for this projection does not
establish arbitrary-polarization accuracy or model a real objective's vector
polarization transfer.

## Discrete operator and three algorithms

`odt_operator` evaluates the centered, padded unitary Fourier transform with
trilinear interpolation at the original q coordinates. Its adjoint is the
exact complex transpose, followed by taking the real part for a real χ.
The separable implementation omits only frequencies never interpolated;
it is tested against the full 3-D FFT, not a reduced physical model.
Illumination weights approximate uniform pupil-area integration using the
forward study's annuli; detector frequency samples have a common pitch.

Write Aχ = c√W χ̂_h and y = c√W g. Internal normalization bounds the operator
norm; `sample_scale` converts back to physical Fourier integrals.
Padding improves interpolation without changing detector data or the native
direct-reconstruction grid.

1. **Direct:** trilinear adjoint deposition on the native volume frequency
   grid, density normalization, inferred Hermitian completion, zero fill,
   inverse FFT. This is an approximate baseline, not A's pseudoinverse.
2. **GP:** a zero-start CGLS minimum-norm least-squares data correction followed
   by χ ≥ 0. Unlike Cartesian measured-bin replacement, this uses the actual
   nonuniform operator. The CGLS normal-residual criterion is not a bound on
   projection error in an ill-conditioned system. For inconsistent data,
   alternating this least-squares projection with positivity is not equivalent
   to minimizing nonnegative least squares; its data residual need not decrease.
3. **Lim TV:** χₖ₊₁ minimizes α TV(χ) + ½‖Aχ − gₖ‖² subject to χ ≥ 0;
   then gₖ₊₁ = gₖ + y − Aχₖ₊₁. Gradient and positivity splits follow Lim
   Appendix A, using isotropic shrinkage and PCG for the actual A* A system.
   `outer_iter=1` gives fixed-data penalized TV, not residual reinjection.

TV(χ) = Σⱼ √[(Dₓχ)ⱼ² + (Dᵧχ)ⱼ² + (D_zχ)ⱼ²], with periodic forward
differences in voxel units. The physical approximation is dx² TV(χ).
The experiment sets α = λ c² dx², λ = alpha_relative × max|χ_direct|;
this preserves the physical data/TV balance when padding or dx changes.
The multiplier is a declared regularization choice, not a fitted true-shape
parameter. Splitting penalties affect convergence, not the objective.
The experiment's `rho` sets both splitting penalties; comparing this parameter
also tests sensitivity to the finite numerical stopping tolerances.

Lim's method is reproduced at the equation level, not by copying its fixed
iteration counts. The native grid, nonuniform data and PCG solve differ from
a diagonal Cartesian-mask implementation. No known spheroid support,
homogeneity, upper RI bound or learned prior is imposed; the finite
computational box and positivity remain assumptions.

## Diagnostics and use with other data

`info.converged` is a numerical criterion, not an accuracy certificate.
TV checks normalized primal/dual residuals, PCG success and stationarity at
the returned nonnegative variable. Outer stopping also requires the declared
data discrepancy unless only one outer solve was requested.
GP checks iterate change and its approximate projection criterion.
Iteration caps remain explicitly unconverged; data fit is reported separately.

```matlab
data = odt_prepare(U,u0,geom,'rytov'); % U is the TOTAL measured field
[n,chi,info] = odt_reconstruct(data,'tv', ...
    struct('alpha',1e-6,'max_iter',200,'outer_iter',5));
```

Use `'born'` for Born data. Geometry can come from a forward call on an
all-water volume; preparation never reads a true spectrum or true phase.
Set `geom.detector_dx` when detector and volume pitches differ.
`demo_inverse(1,'rytov')` provides a same-model voxel demonstration, not the
independent reference experiment; its default parameters are illustrative.

RI error is ‖n_rec − n_true‖₂ / ‖n_true − nₘ‖₂ over the entire volume,
not divided by the large background RI. Physical-field error compares
predicted and exact co-polarized **scattered** fields with common angle weights.
Axial span uses the true half-contrast threshold on the central line;
the voxelized truth spans 5.9 μm, versus the analytic diameter 6 μm.
`axial_span_censored` marks profiles that reach a volume boundary; their span
is window-limited, not a resolved object thickness.
Small field error does not ensure small RI error in a missing-cone problem.

## References

- Lim et al. (2015), *Comparative study of iterative reconstruction algorithms
  for missing cone problems in optical diffraction tomography*, Eqs. 14–15
  and Appendix A, [local PDF](../ref/background/oe-23-13-16933.pdf),
  DOI 10.1364/OE.23.016933.
- Goldstein and Osher (2009), *The Split Bregman Method for L1-Regularized
  Problems*, [local PDF](../ref/background/goldstein2009.pdf),
  DOI 10.1137/080725891.

Historical low-NA/nearest-bin results in `verification.md` and provisional
runs in `results/provisional/` are not the current experiment's conclusions.
