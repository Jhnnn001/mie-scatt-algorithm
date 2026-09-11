# mie-scatt-algorithm 프로젝트 보고서

작성일: 2026-09-13
범위: 2026-09-11 ~ 2026-09-13에 이 저장소에서 진행한 (1) spheroid 정확해 solver, (2) Born/Rytov forward 모델과 오차 분석, (3) ODT inverse 복원. 각 폴더의 `plan.md`, `history.md`, `readme.md`, `report.md`, CSV, log를 근거로 정리했다. 수치는 모두 저장된 파일에서 가져온 값이며, 이 보고서를 위해 새로 계산한 값은 없다.

표기: 수식은 LaTeX 대신 일반 문자와 기호로 쓴다. 여기서 "exact"는 닫힌 형태의 해가 아니라, 모드 절단·경계 잔차·cutoff 검사를 통과한 수치적 full-wave reference를 뜻한다. Scalar exact(scalar Helmholtz 경계값 해)와 Maxwell exact(full-vector 해)는 서로 다른 방정식의 해다.

---

## 0. 목적

최종 목표는 optical diffraction tomography(ODT) 스타일의 forward/inverse 모델링이다. 실제 실험(Tomocube 계열 Rytov 기반 ODT, 532 nm, 물 배경, 세포 굴절률 약 1.37)과 같은 조건에서 다음 질문에 답하려 했다.

1. 세포 크기의 균질 spheroid에 대해 Maxwell 방정식의 정확한 산란장을 계산할 수 있는가. 이것이 모든 근사법의 ground truth가 된다.
2. 실제 ODT가 쓰는 scalar Born/Rytov 근사가 이 ground truth에 대해 얼마나 정확한가. 특히 논문 굴절률(1.37)과 세포 크기 구조에서 Rytov가 왜 깨지는가.
3. Rytov 데이터를 이용한 굴절률 복원(direct Fourier, Gerchberg–Papoulis, TV)이 실제 산란장 데이터에서 어떤 결과를 주는가. missing cone과 근사 오차가 복원에 어떻게 나타나는가.

세 축은 순서대로 의존한다. exact solver가 없으면 forward 근사의 오차를 정량화할 수 없고, forward 오차를 모르면 inverse 결과에서 수치 오차·근사 오차·missing cone 효과를 분리할 수 없다.

---

## 1. 전체 흐름 요약

| 시각 (KST) | 작업 | 상태 |
|---|---|---|
| 09-11 23:10 | 참고 논문(Asano & Yamamoto 1975, Barton 2001) 추가, 폴더 생성 | commit `bacc44d`, `97618ee` |
| 09-12 04:36 | spheroid solver plan v1, 좌표·평면파·Gauss–Legendre 구현 | commit `5cecb3f` |
| 09-12 04:12–13:53 | spheroidal special functions, outgoing ODE, Galerkin solver, 10개 test suite 통과 | commit `569307b` |
| 09-12 15:38–17:49 | biological baseline(prolate 5×2.5 μm, n=1.37)에서 scalar-polarization 특성화 | commit `b0d937b`, `51094f1` |
| 09-12 | forward plan, `born_rytov_spheroid.m`(surface-integral), `born_rytov_voxel.m`(FFT/Ewald) 구현 | untracked `forward/` |
| 09-12 | inverse v1 (linear ODT, nearest-bin gridding, direct/GP/TV) 구현·검증 | untracked `inverse/` |
| 09-12 18:45–18:56 | voxel Rytov의 detector-NA 적용 순서 버그 발견·수정, 15-case 비교 | `forward/tests/rytov_analysis.*` |
| 09-12 19:10–19:24 | scalar prolate 경계값 reference로 scalar/Rytov 오차 분리, Q 항·field zero·geometry control | `forward/exp/report.md` |
| 09-12 | oblate(3×5 μm, n=1.365) 구조로 재설계, scalar-vs-vector NA pilot, 81-direction Born/Rytov sweep, xz 그림, Q 메커니즘 | `forward/exp/oblate_*`, `forward/history.md` |
| 09-13 | 크기 실험(폴더 `1`), 굴절률 대비 실험(폴더 `2`), Δn/100 xz 그림 | `spheroid-analytic-forward/1`, `2` |
| 09-13 01:11–01:43 | inverse v2: matched Ewald operator, exact Maxwell 데이터로 direct·GP 복원 완료, TV·padding 8 실행 중 | `inverse/results/` |

핵심 결론을 먼저 요약하면 다음과 같다.

- **Exact solver**: 균질 prolate/oblate spheroid의 full-vector 산란장을 계산하는 MATLAB solver가 완성되었고 Barton 2001 benchmark, 광학정리, 흡수 에너지 수지, 경계 연속성 등을 통과했다.
- **논문 굴절률 + prolate 구조에서 Rytov가 깨진 이유**: n=1.37, z 방향 10 μm 두께에서 1차 Rytov가 버리는 항 Q = ∇ψ·∇ψ(로그장 기울기의 제곱)가 작지 않다. 중심선 위상 지연 4.09 rad 자체보다, 중앙(두꺼움)과 가장자리(얇음)의 위상 차이가 좁은 횡방향 거리(b=2.5 μm) 안에서 생기면서 만드는 큰 위상 기울기, 내부 반사파가 축 근처에 모여 만드는 급격한 위상 변화, 그리고 검출면 바로 앞의 field zero가 원인이다. Scalar 근사 자체의 오차는 0.25%에 불과했다.
- **NA와 FFT**: Rytov에서 detector pupil은 지수화 **뒤에** 적용해야 한다(74% → 26%). Scalar 근사의 정확도는 detector NA와 입사각에 따라 달라지며, 작은 NA가 항상 유리하지 않다. voxel FFT 모델의 Ewald 보간/padding 오차는 수 % 수준으로 무시할 수 없어 연속 형상의 해석적 Fourier transform을 대조군으로 두었다.
- **구조 재설계**: oblate 3×5 μm, n=1.365로 바꾸면 Rytov 오차가 25.7% → 10.5%로 줄지만 5% 목표에는 못 미친다. Δn을 1/100로 줄이면(n_p=1.33568) scalar exact 대비 0.09%, Maxwell co-pol 대비 0.66%로 떨어진다. 이 조건이 inverse 실험의 입력이 되었다.
- **Inverse**: matched operator를 이용한 direct 복원은 굴절률 대비 상대 오차 약 47%, 축 방향 두께 4.3 μm(참값 5.9 μm)로 missing cone 효과가 뚜렷하다. 비음수 GP는 36.5%로 줄이지만 100 iteration에서 미수렴이다. TV와 padding 8 대조 실행은 이 보고서 작성 시점(09-13 01:43)에 진행 중이다.

---

## 2. 축 1: spheroid의 exact 계산 (`spheroid-analytic-forward/`)

### 2.1 목표와 방법

목표는 `[Ex,Ey,Ez,info] = spheroid_field(lambda,n_p,n_m,a,b,theta_inc,phi_inc,pol,X,Y,Z,opts)` 하나로 임의의 입사 방향·복소 편광·관측점(내부/외부)에서 복소 전기장을 얻는 것이다. `a`는 z 반축, `b`는 xy 반축이며 a>b prolate, a<b oblate, a=b는 Mie로 분기한다. 시간 convention은 exp(−iωt), 흡수(복소 n_p)를 허용한다.

방법은 spheroidal coordinate separation of variables다.

- 스칼라 파동함수 Π_lm = S_lμ(c;η) R_lμ(c;ξ) e^{imφ}, 벡터 파동함수 M = ∇Π × r, N = (1/c)∇×M.
- 외부 산란장은 outgoing radial 함수, 내부장은 regular radial 함수로 전개하고, 표면에서 접선 E·H 연속 조건을 Barton 2001 방식의 Galerkin(표면 투영) 으로 맞춘다. 입사장은 표면 노드에서 수치 투영한다.
- 각 azimuthal order m마다 4(L−μ+1)개 미지수의 선형계를 푼다.

Asano & Yamamoto 1975는 해석적 경계 매칭과 far-field 식(Mie의 π_n, τ_n의 spheroidal 대응)을, Barton 2001은 임의 입사장에 대한 Galerkin 투영과 near-field 식을 제공한다. 두 논문 모두 `ref/mie-scatt/`에 있고, spheroidal 함수의 recurrence는 DLMF Ch. 30과 Flammer 1957을 따랐다.

### 2.2 수치적으로 어려웠던 점과 해결

공식 자체보다 부동소수점에서 안정하게 만드는 일이 훨씬 어려웠다. 세 차례 plan 검토와 여러 hardening 단계에서 다음이 확정되었다.

**좌표 변환 (`sph_coords.m`)**
- 축과 초점 근처에서 ξ²−1, 1−η²를 큰 수의 차로 계산하면 유효숫자를 잃는다. α=ξ²∓1, β=1−η², D=ξ²∓η²를 상쇄 없는 식으로 직접 계산해 반환하고, 후속 함수는 이를 재계산하지 않는다.
- 초점(prolate)과 focal ring(oblate)은 좌표계의 singular set이다. 거리 δ=1e−8 안의 점은 방위각과 접근 side를 보존한 채 δ 밖으로 투영한다. 축, focal segment, focal disk는 regular limit로 처리하고 옮기지 않는다.
- oblate z=0에서 sign(0):=+1로 sheet를 고정하고, 반대 sheet를 골라도 parity 때문에 같은 장이 나옴을 확인했다.
- 나중에 oblate xz 그림 계산 중 축 위 네 점에서 |η|가 반올림으로 1을 넘는 문제가 발견되어 `eta_abs = min(1, eta_abs)` 한 줄을 추가했다(uncommitted 변경, regression test 추가).

**Outgoing radial 함수 (`sph_radial_ode.m`)**
- 초기 설계의 spherical Neumann series R^(2)는 항 비율이 1/ξ²이라 ξ≈1이나 oblate ξ<1에서 발산·overflow했다.
- 최종 설계는 canonical normalization이 경계 매칭에 불필요하다는 점을 이용한다. 각 mode마다 R_out(ξ0)=1로 임의 정규화하고, inverse-ξ asymptotic seed(최대 80항)에서 Riccati/log-amplitude ODE(`ode113`, RelTol 1e−12)를 안쪽으로 적분한다. Galerkin 계수가 임의 scale을 흡수한다.
- 두 독립 시작점의 log-derivative·field 일치(<1e−10)와 ODE/asymptotic backward error(≤1e−9)를 검사하고 실패하면 해당 basis를 거부한다(fail-closed).

**Angular 함수 (`sph_angular.m`)**
- 정규직교 Ferrers basis에서 도함수는 degree-dependent square-root factor가 필요하다. raw S′, S″ 대신 벡터장이 실제로 요구하는 metric-weighted 조합 W0=S, W1=√β S′, W2=β^{3/2} S″, Wm=μS/√β를 직접 계산해 축에서 0/0을 없앴다.

**Eigenproblem·특수함수**
- prolate c=50의 일부 regular radial mode는 normalization sum의 조건수가 1e15에 달해 Frobenius initial-value fallback을 도입했다. 표현 민감성이 남는 high-order complex oblate(예: c=200+20i, L=234)는 조용히 틀린 값을 내지 않고 error를 낸다.
- MATLAB `besselj`/`besselh`가 order×argument 외적 broadcasting을 하지 않아 배열을 명시적으로 구성했다.

**Galerkin 선형계 (`spheroid_solve.m`)**
- 초기 rcond가 1e−20까지 떨어졌다. column equilibration으로 완화하고, 축방향 입사에서 해석적으로 존재하지 않는 sector만 정확히 0으로 두며, rank-deficient면 `lsqminnorm`을 쓴다.
- 수렴 루프: incident spectrum tail, Q_φ 독립성, 상위 세 l shell의 field 기여, 독립 grid의 경계 잔차, Q_θ 재계산, (L+5, M+2) 재계산, 특수함수 진단을 모두 통과해야 `info.validated=true`.

**Runtime-validation envelope**: x_ext=k_m·max(a,b) 0.02–250, x_int ≤250, a/b 1/4–4, Re(n_p/n_m) 0.5–2, Im(n_p)/Re(n_p) ≤0.1, M≤300, L≤400. 이 범위는 자격 조건이며, 매 실행의 검사 통과가 별도로 필요하다.

### 2.3 검증

commit `569307b`에서 10개 suite 통과(로그: `../tmp/spheroid-analytic-forward/green/run_all_final10.log`).

| 검사 | 결과 |
|---|---|
| Wiscombe MVTstNew case 14 Mie benchmark, Rayleigh limit, 광학정리 | 통과 |
| Barton 2001 Fig. 7/13 far-field sampled maxima (prolate·oblate 2:1, n=1.33/1.50) | 상대오차 7.92e−8, 5.66e−7, 3.02e−7, 4.25e−7 (기준 1e−4) |
| lossless 에너지 수지 C_ext=C_sca | prolate 9.89e−8, oblate 4.43e−8 |
| lossy: Poynting surface flux와 volume loss로 각각 구한 C_abs | 1e−5 안에서 일치 |
| 경계 연속성(접선 E, H, n²E_ξ) | 최대 1.86e−7 |
| curl M=cN, curl N=cM, div=0 유한차분 | 통과 (축·focal segment·disk 포함) |
| 길이 1e±200배, 진폭 1e±200배 불변성 | 통과 |

실행시간: `test_physics` 7,520.8 s, `test_barton2001` 1,594.5 s. 537-case stress campaign은 정의와 resume 메커니즘만 smoke test했고 전체 수치 실행은 하지 않았다.

**Scalar-polarization 특성화 (T10, `test_scalar_approximation.m`)**: 세포 baseline(λ=0.532 μm, n_m=1.335381534, n_p=1.37, prolate 5×2.5 / oblate 2.5×5 μm, x_ext=78.86)에서 산란장을 입사 편광 e0에 투영했을 때 버려지는 성분의 비율 ε_perp를 쟀다. 좁은 forward patch(|u|,|v|≤0.08)에서 0.020–0.035, 넓은 patch(≤1)에서 0.126–0.203. 6 case에 2시간 11분 소요. 좁은 각도에서만 fixed-polarization scalar 축약이 성립한다는 뜻이며, 넓은 각도 leakage에는 far-field transversality(E_sca·r̂=0)가 강제하는 기하학적 투영도 들어 있다.

### 2.4 Scalar 경계값 reference (실험용, `forward/exp/scalar_prolate.m` 및 batch)

Born/Rytov 오차를 "scalar 근사 오차"와 "Rytov 절단 오차"로 분리하려면 scalar Helmholtz 방정식의 exact 해가 따로 필요하다. 같은 spheroidal special function(eigen, angular, radial, outgoing ODE)을 재사용하되, Maxwell 경계조건 대신 U와 ∂ₙU의 연속성을 맞추는 scalar 경계값 solver를 조립했다.

- 처음 버전은 real-index prolate, +z 축방향 입사, m=0에 한정.
- 이후 oblate·tilted incidence용 batch reference(`forward/exp/oblate_na_sweep.m`)로 확장: 각 θ에서 xz 평면 TM/TE 두 RHS를 풀고, signed m마다 Galerkin matrix를 공유해 여러 입사각을 한 번에 푼다. coarse (L,M)=(109,55)와 fine (114,59)를 비교한다.
- 검증: near-sphere에서 independent scalar spherical partial-wave 해와 1.12e−6, 경계 잔차 5.47e−15, L=112 vs 132 far amplitude 1.21e−12, weak-potential 미분과 analytic Born 계수 5.10e−5.

이 reference는 public API가 아니라 원인 분리용 실험 코드이며, vector solver와 특수함수를 공유하므로 완전히 독립된 구현은 아니다.

### 2.5 남은 한계

- 독립 문헌 비교는 Barton의 네 far-field intensity maximum뿐이다. 복소 위상·near field에 대한 독립 spheroid benchmark는 없다.
- 평면파 입사만 지원한다. focused beam은 입사장 투영만 바꾸면 되지만 구현하지 않았다.
- ξ0=0인 극한 oblate(원판)는 outgoing log-ODE normalization 때문에 지원하지 않는다.
- 11개 suite(T10 포함)를 한 tree에서 한 번에 통과시킨 로그는 아직 없다.

---

## 3. 축 2: forward — Born/Rytov 근사와 오차 분석 (`forward/`)

### 3.1 두 가지 forward 구현

**(a) 연속 형상 surface-integral 모델 (`born_rytov_spheroid.m`)**: 균질 spheroid의 Born 체적 적분

u_B(r) = f ∫_V G(|r−r′|) u0(r′) dV′,  G(R)=e^{ikR}/(4πR),  f = k0²(n_p²−n_m²)

을 Green 제2 항등식과 보조함수 w=(k̂·r′)u0/(2ik)로 **표면 적분 + jump 항**으로 정확히 바꾼다.

u_B(r) = f { ∮_S [G ∂ₙw − w ∂ₙG] dS′ − χ_V(r) w(r) }

내부 관측점의 1/R singularity가 사라지고 Gauss–Legendre(θ)×trapezoid(φ)로 spectral 수렴한다. 먼 점에서는 공통 위상 e^{ikρ}를 분리하고 R−ρ=(|r′|²−2r·r′)/(R+ρ)로 cancellation-free하게 계산한다(|r|=1e9a에서 naive kernel 오차 2e−8 → 3e−15). Rytov는 u_R = u0·expm1(u_B/u0). tol 기준 coarse/fine quadrature 쌍으로 Born·Rytov를 각각 검사한다. Sphere partial-wave series, brute-force volume quadrature, Rayleigh–Gans far field에 대한 test(T1–T5)를 구현 당시 통과했다. 다만 2026-09-13의 별도 실행에서 `forward/tests/run_all.m`은 `test_born_rytov.m:178`에서 멈춘다. test는 `tol=99*eps`에 대해 `InvalidTolerance` error를 기대하지만 현재 함수는 양수 tol만 검사한다. 이 함수/test 불일치는 아직 고치지 않았다. 이후 실험은 모두 `born_rytov_voxel.m` 또는 실험용 batch reference를 사용하며 이 legacy 함수를 호출하지 않는다.

**(b) Voxel FFT/Ewald 모델 (`born_rytov_voxel.m`)**: 3차원 굴절률 배열 n(x,y,z)를 받아 +z detector 평면의 Born/Rytov 산란장을 계산한다. 실제 ODT forward의 구조와 같다.

1. f = k0²(n²−n_m²)를 zero-padding하고 3D FFT × dx³.
2. 각 입사 방향에 대해 q = k_s − k_i인 shifted Ewald hemisphere에서 spectrum을 linear 보간.
3. propagating hemisphere(kz>0)에 i·exp(i kz z_det)/(2kz)를 곱한다.
4. 2D inverse FFT로 detector 앞의 linear Born field를 만든다.
5. Rytov: ψ₁=B_pre/U0, U_R = P[U0·expm1(ψ₁)], Born: U_B = P[B_pre]. P는 detector pupil.

Evanescent와 정확히 grazing인 mode는 제외한다. 입력 검사로 dx < λ/(2n_m) (=0.199194 μm), 입사 이동된 Ewald sample의 FFT 대역폭, 입사 reference가 pupil 안에 있는지, z_det가 voxel grid +z 밖인지를 확인한다.

두 구현 모두 scalar이며 polarization 인자가 없다. Vector reference와 비교할 때는 co-polarized 투영 u_co = E_s·conj(e0)를 쓴다.

### 3.2 논문 굴절률 + prolate baseline에서 Rytov가 깨진 과정

Tomocube 계열 논문(`ref/tomocube-rytov/`, Kim et al. 2014)의 조건을 그대로 가져왔다.

| 항목 | 값 |
|---|---|
| λ₀ / n_m / n_p | 0.532 μm / 1.335381534 (20 °C 물) / 1.37 |
| 형상 | prolate a=5 μm (z), b=2.5 μm. 전체 5×5×10 μm, 긴 축이 입사 방향 |
| 입사 / 편광 | +z / x |
| detector | z=7 μm (끝에서 2 μm), NA=0.1, 128³ voxel, dx=0.1 μm, FOV 12.8 μm |
| 중심선 ray 위상 지연 k0·Δn·2a | 4.0886 rad |

Detector pupil 안 복소 산란장의 상대 L2 오차 ‖model−exact‖/‖exact‖ (exact = validated Maxwell 해의 Ex 투영):

| 모델 | 오차 |
|---|---:|
| Born | 197.56% |
| 수정 전 voxel Rytov | 74.40% |
| 수정 전 Rytov에 final pupil만 추가 | 26.21% |
| **수정된 voxel Rytov** | **25.74%** |
| 연속 형상 Rytov (같은 detector grid) | 24.80% |
| 연속 형상 Rytov, window 51.2–102.4 μm | 24.30–24.34% |

즉 Born은 완전히 실패하고, 올바르게 구현한 1차 Rytov도 약 25%를 남긴다. 상대 굴절률 차 0.001 대조군에서는 1.21%(voxel)/1.04%(연속)로 떨어지므로 대비에 민감하다. voxel/연속 field 차이는 padding 2에서 3.41%, padding 4에서 1.44%이므로 voxel화나 유한 window로 24%를 설명할 수 없다.

**Scalar 근사 vs Rytov 절단의 분리** (scalar prolate 경계값 reference, NA=0.1, z=7 μm):

| 비교 | 오차 |
|---|---:|
| Scalar exact vs Maxwell Ex | 0.2527% |
| 1차 Rytov vs scalar exact | 24.174% |
| 1차 Rytov vs Maxwell Ex | 24.302% |

따라서 이 조건에서 25%의 대부분은 **scalar 축약이 아니라 1차 Rytov 절단** 때문이다.

### 3.3 왜 깨졌나: 생략된 항 Q와 "중앙의 phase 차이"

로그장 ψ = ln(U/U0), U0 = exp(i k_i·r)로 두면 scalar Helmholtz 방정식은 정확히

∇²ψ + 2i k_i·∇ψ + Q = −f,  Q = ∇ψ·∇ψ = (∂ₓψ)² + (∂ᵧψ)² + (∂_zψ)²

이다. Q에는 복소켤레가 없다. 1차 Rytov는 Q를 버리고 ψ₁ = u_B/u0로 둔다. ∇²ψ는 남기므로 회절을 버리는 근사가 아니다. ψ=α+iΔΦ(로그 진폭 α, 상대 위상 ΔΦ)로 쓰면

Q = |∇α|² − |∇ΔΦ|² + 2i(∇α·∇ΔΦ)

이므로 Q가 작으려면 **위상의 공간 기울기**와 진폭 기울기가 모두 작아야 한다. 큰 누적 위상을 지수로 표현하는 것과, 그 위상이 공간적으로 급격히 변하는 것을 정확히 표현하는 것은 다른 문제다. Q와 f의 단위는 모두 길이⁻²이므로 국소 비교 기준은 |Q| ≪ |f|이며, 이번 물체에서 |f|≈11.16 μm⁻²(oblate) 또는 그 이상이다.

**왜 중앙에서 위상 차이가 커지는가 (ray 직관)**: 두께 t(ρ)=2a√(1−ρ²/b²)인 spheroid를 지나는 ray의 위상은 φ(ρ) ≈ k0·Δn·2a√(1−ρ²/b²)다. 중앙에서 4.09 rad, 가장자리에서 0이며 이 차이가 b=2.5 μm 안에서 생긴다. ρ=b/√2에서 기울기는 1.635 rad/μm, |∇φ|²/f ≈ 0.205다. 두께가 연속이어도 그 기울기는 bounded되지 않는다. 즉 "중심선 위상이 4 rad"가 아니라 "중심과 가장자리의 위상 차이가 좁은 횡방향 거리에서 생긴다"가 문제다.

이 직관을 검증하는 **geometry control** (같은 4.0886 rad 위상, 같은 z=7 μm, NA=0.1):

| 형상 | 1차 Rytov 산란장 오차 |
|---|---:|
| prolate b=2.5 μm | 24.17% |
| sphere b=5 μm | 12.72% |
| 무한 평판(두께 10 μm) | 2.98% (전체장 5.30%) |

같은 중심 위상이라도 횡방향 구조가 넓을수록 오차가 작다. 중심선 위상 하나로 Rytov 정확도가 결정되지 않는다.

**Q를 직접 계산한 결과** (scalar exact의 해석적 미분 ∇ψ=∇U/U−ik_i 사용, branch/unwrap 무관):

| 영역 (prolate, t=1) | 체적 가중 RMS \|Q\|/f, dx 0.05 μm | dx 0.025 μm |
|---|---:|---:|
| 중앙 원통 ρ<b/2 | 0.2410 | 0.2394 |
| 바깥 ρ>0.8b | 0.2863 | 0.2858 |

potential을 1/10로 줄이면 0.0149, 0.0330으로 떨어진다. 횡방향 성분이 norm의 91–99%를 차지한다.

**Field zero**: scalar exact의 pupil 이전 전체장에 ρ=0.1598748 μm, z=6.9963225 μm에서 zero가 있다(|U|=2.41e−14 at L=112, 8.64e−12 at L=132). 축대칭이므로 ring zero이며 검출면 z=7 μm 바로 0.0037 μm 앞이다. 그곳에서 ψ=ln(U/U0)는 singular하므로 "작은 log-gradient" 가정이 성립할 수 없다.

**검출 거리 의존성과 전파 비일관성**: 1차 Rytov를 각 검출면에서 새로 평가하면 z=5.2/7/14 μm에서 16.87%/24.17%/47.21%다. z=5.2 μm에서 만든 Rytov 장을 물속에서 선형 전파한 것과의 차이는 z=7에서 6.53%, z=14에서 26.93%. 전파된 선형 위상을 지수화하는 것과 지수화된 장을 전파하는 것은 다르며, 이 차이가 곧 생략된 Q 항의 효과다.

**2차 로그항 진단**: 약한 ±h potential에서 exact scalar 해를 대칭 차분해 ψ₂ = U₂/U0 − ψ₁²/2를 독립 추출하고 U0·exp(tψ₁+t²ψ₂)를 평가했다. z=7 μm에서 24.17% → 13.38%, z=5.2 μm에서 16.87% → 4.39%. 그러나 z=14 μm에서는 2,705%로 발산하므로 2차항이 일반적 해결책은 아니다.

**대비 의존성** (z=7 μm, potential scale t):

| t | Born | 1차 Rytov | 2차 진단 |
|---:|---:|---:|---:|
| 0.10 | 16.49% | 2.52% | 0.046% |
| 0.25 | 41.37% | 6.15% | 0.265% |
| 0.50 | 84.86% | 12.01% | 1.233% |
| 1.00 | 197.23% | 24.17% | 13.381% |

### 3.4 Detector NA 적용 순서 버그와 수정

수정 전 `born_rytov_voxel.m`은 detector NA로 잘라낸 Born 장으로 Rytov 위상을 만들었고, 지수화 뒤 pupil 밖에 다시 생긴 주파수도 출력에 남겼다. 선형 pupil과 비선형 지수화는 교환되지 않는다. 물리적으로 1차 Rytov 위상은 objective 앞의 unpupilled propagating field에서 만들고, 지수화 후 마지막에 pupil을 적용해야 한다.

수정 효과: 74.40% → 25.74%. 이 중 대부분은 pupil 밖 주파수 제거(26.21%)이고 순수한 순서 수정은 26.21% → 25.74%다. FOV 25.6 μm에서는 final-pupil-only 41.98% vs 올바른 순서 25.33%로 차이가 커진다. 회귀 test(`test_born_rytov_voxel.m`)에서 detector NA 변경 시 pupil 안 spectrum 변화 2.1e−16, pupil 밖 누설 1.9e−16. 부작용으로 `info.rytov_phase`는 pupil 이전 ψ₁이며 pupil 이후 log(1+U_R/U0)와 다르다. 이 수정 때문에 기존 inverse v1 test(dx=0.2)는 `RytovNyquistViolation`으로 멈추게 되었다.

### 3.5 Scalar 근사와 NA: 얼마나 넓은 NA까지 scalar를 써도 되는가

Oblate 구조(3.6절)에서 scalar exact와 Maxwell exact를 detector NA에 따라 비교했다(console 기록, `spheroid-analytic-forward/history.md` §13). 정상입사, 모든 값 %:

| Detector NA | ε_perp (편광 밖) | ε_co (co-pol 진폭·위상 오차) | ε_full (전체 벡터) |
|---:|---:|---:|---:|
| 0.10 | 1.902 | 0.076 | 1.904 |
| 0.50 | 2.715 | 0.303 | 2.732 |
| 1.00 | 3.024 | 0.726 | 3.110 |
| 1.30 | 3.263 | 1.468 | 3.578 |

전체 벡터 오차의 1% crossing은 NA≈0.0281, co-pol만의 1% crossing은 NA≈1.1755다. "scalar 근사가 성립한다"는 말에는 어떤 성분과 허용 오차를 말하는지 반드시 붙어야 한다.

Illumination NA 0–0.5(12개 극각, φ=0/45/90°)로 확장하면 NA_det=0.1에서 최악 ε_full이 37.99%, NA_det≥0.5에서는 801개 표본 최대 3.697%다. **작은 detector NA가 비정상입사에서는 오히려 불리하다**: 기울어진 강한 forward lobe가 고정된 +z pupil을 벗어나므로 좁은 pupil은 약한 큰 각 산란만 모아 상대 오차가 커진다. NA_det<NA_illum이면 직접 입사광도 pupil 밖이다. 이 결과가 detector NA=1.0을 최종 실험 조건으로 정한 근거다.

### 3.6 구조 재설계: oblate 3×5 μm, n=1.365 (81-direction sweep)

Prolate baseline의 25%를 줄이기 위해 z 방향으로 납작한 oblate(핵 부피 약 300 μm³ 수준)로 바꾸고 대비를 약간 낮췄다.

| 항목 | prolate baseline | 새 oblate |
|---|---:|---:|
| n_p | 1.370 | 1.365 |
| a (z), b (xy) | 5, 2.5 μm | 3, 5 μm |
| 전체 x×y×z | 5×5×10 | 10×10×6 μm |
| detector z / NA | 7 μm / 0.1 | 5 μm / 1.0 |
| illumination NA | 0 | 0–0.5 (θ≤21.99°), 81 방향 |
| 중심선 위상 | 4.0886 rad | 2.0989 rad |

편광은 기존 `pol=[1,0]` convention e0 = cosφ·e_TM − sinφ·e_TE = Rz(φ)Ry(θ)Rz(−φ)x̂를 유지했다. 축대칭 물체라도 이 정의는 φ에 따라 TE/TM 혼합비가 달라지므로 각 θ에서 TM·TE를 모두 풀고 회전·선형결합으로 복원했다. 입사각 평균은 illumination pupil **면적 균일 가중**(annulus 면적/방위각 수)이다.

결과 (co-polarized Maxwell reference, pupil 이후 복소 산란장, %):

| 모델 | 가중 평균 | pooled L2 | 최소 | 최대 |
|---|---:|---:|---:|---:|
| 연속 형상 Born | 91.319 | 91.368 | 89.084 | 93.648 |
| **연속 형상 Rytov** | **10.509** | 10.528 | 9.869 | 11.192 |
| voxel Born (256×256×80, pad 2) | 89.679 | 89.716 | 88.095 | 91.968 |
| voxel Rytov | 9.751 | 9.775 | 8.500 | 10.844 |
| Scalar exact (scalar 근사 자체) | 0.705 | 0.705 | 0.658 | 0.733 |

- Scalar exact 기준으로도 Born 91.281% / Rytov 10.525%이므로 잔차는 scalar 축약이 아니라 Rytov 절단이다.
- Rytov는 25.7% → 10.5%로 개선되지만 5% 목표에는 못 미친다.
- illumination NA에 따라 Rytov 오차가 9.869%(0) → 11.187%(0.5)로 완만히 증가한다. 같은 NA에서 방위각 차이는 최대 0.0091 pp.
- voxel의 9.751%가 연속 형상보다 작은 것은 수치 오차와 근사 오차의 부분 상쇄다. voxel-연속 field 차이는 padding 2에서 4.505%, padding 3에서 2.027%로 보간 오차가 무시할 수 없다.
- 수치 검증: cutoff (109,55) vs (114,59) 6.5e−13, 독립 경계 잔차 4.7e−12(vector)/1.9e−12(scalar), window 25.6→51.2 μm 변화 ≤0.0049 pp, dx 0.1→0.05 변화 0.0001%.

### 3.7 xz 그림과 지표 혼동(10% → 5%)의 정정

θ=0/5/10/15°, y=0 단면에서 전체 |Ex|/|E0|(입사장 포함, pupil 없음)를 Exact/Born/Rytov 4×3 그림으로 저장했다(`forward/exp/oblate_Ex_xz_4x3.png`). 정상입사 최대 |Ex|는 Exact 1.470, Born 2.774, Rytov 1.590으로 Born이 downstream 진폭을 크게 과대평가한다.

같은 xz 배열의 오차: Rytov 복소장 8.3–9.0%, |Ex| 크기만 5.2–5.8%, scalar vs Maxwell 복소 2.24–2.69%. "10%가 5%로 줄었다"는 인상은 알고리즘 개선이 아니라 (i) pupil 이후 산란장 → raw 전체장, (ii) 복소 → 크기만, (iii) 81 방향 평균 → 4 각도라는 **지표 변경**이었다.

### 3.8 Q는 어디에서 커지는가: 꼭지점이 아니라 내부 중앙

Oblate exact scalar의 내부 xz 단면에서 |Q|/f 통계(θ=0°): median 0.038, 90th percentile 0.197, 99th 1.74, RMS 0.478, >1인 점 1.686%. 대부분 작지만 좁은 띠에서 매우 크다. Sampled maximum은 θ=0°에서 (x,z)=(0,−0.20) μm, |Q|/f=14.33이며, z=±3 μm 꼭지점이 아니라 **내부 중앙**이다. θ가 커지면 띠가 +x로 이동한다. |x|<0.5, |z|<2 μm 영역이 Σ|Q|²의 95.18%를 차지한다.

최대점의 항 분해(f로 나눔): |∇α|² = +2.82, −|∇ΔΦ|² = −11.50, 2i∇α·∇ΔΦ = +11.40i. z 방향 상대 위상 기울기 11.33 rad/μm, 로그 진폭 기울기 5.61 /μm. |U|=0.704로 0에 가깝지 않은데도 Q가 크다.

**메커니즘**: 상부 곡면이 내부에서 보면 오목 반사면 역할을 하여 내부 반사파가 축 근처로 모인다(caustic). 증거는 세 가지다.
- x=0에서 Hann-windowed z-FFT: +kz peak 16.12 μm⁻¹(≈k_p=16.12), −kz peak −9.23 μm⁻¹, 진폭비 0.110. 음의 kz 성분 비율은 축에서 0.065, x=2 μm에서 0.0003, 순수 평면파 window leakage는 2.2e−8.
- 실제 굴절률로 Snell 굴절 + 1회 내부 반사 ray tracing: 입사 반경 1/2/3/4 μm의 반사 ray가 축과 만나는 z = −1.04, −0.77, −0.24, +0.86 μm. Q 최대 영역과 겹친다.
- 반사파와 전진파의 간섭이 급격한 진폭·위상 변화를 만든다.

Fresnel 가중이나 Debye 반사 차수 분해는 하지 않았으므로, 이 반사가 10% detector 오차의 몇 %를 설명하는지는 배정하지 않는다. Q가 큰 위치와 detector 오차의 기여 위치는 다른 질문이다(Green 전파와 간섭이 개입한다).

### 3.9 2차 Rytov와 대비 sweep (oblate)

Scalar exact 기준, z=5 μm, NA=1, θ=0°: 1차 Rytov 9.890% → 2차 진단 3.847%. potential scale t=0.1/0.25/0.5/1에서 1차 0.867/2.193/4.512/9.890%, 2차 0.034/0.214/0.880/3.847%. 16개 조건 모두 2차가 1차보다 작았다. 고차항 절단이 상당한 원인이라는 판단을 지지하지만, 오차장은 복소 간섭하므로 "60%는 Q 때문"과 같은 가산 분할은 하지 않는다.

### 3.10 크기 실험과 대비 실험: 무엇을 바꿔야 Rytov가 정확해지는가

두 실험은 같은 runner(`spheroid-analytic-forward/1/parameter_sweep.m`)와 81 방향, 같은 오차 정의를 쓴다.

**크기 실험 (`1/`)**: a:b=3:5 유지, 두 반축을 s=1, 0.1, 0.01배. 검출면 z=5 μm 고정.

| s | Scalar Born vs scalar exact | Scalar Rytov vs scalar exact | Vector Born vs Maxwell | Scalar exact vs Maxwell co | Scalar Rytov vs Maxwell co |
|---:|---:|---:|---:|---:|---:|
| 1 | 91.28% | 10.53% | 91.27% | 0.71% | 10.51% |
| 0.1 | 9.43% | 8.90% | 9.22% | 6.72% | 11.40% |
| 0.01 | 0.72% | 0.72% | 0.78% | 30.52% | 29.99% |

크기를 줄이면 Born/Rytov 차수 오차는 줄지만, 작은 물체는 넓은 각도로 산란하므로 **scalar 모델 자체가 깨진다**(0.1 μm 물체에서 30%). 횡방향 투영을 넣은 vector Born 대조군은 0.78%이므로 이 30%는 Born 차수가 아니라 벡터 방향성 누락 때문이다. 산란 신호 norm은 1 → 0.0115 → 3.7e−5로 급감한다.

**대비 실험 (`2/`)**: 형상 유지, Δn을 t=1, 0.5, 0.01배 (n_p = n_m + tΔn, Δn=0.029618466).

| t | n_p | δ (rad) | Rytov vs scalar exact | Rytov vs Maxwell co | Rytov vs Maxwell full | Born vs scalar exact |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.365 | 2.099 | 10.525% | 10.509% | 10.935% | 91.28% |
| 0.5 | 1.350190767 | 1.049 | 4.760% | 4.802% | 5.534% | 44.37% |
| 0.01 | 1.33567771866 | 0.021 | 0.0905% | 0.659% | 2.736% | 0.878% |

t=0.01에서 Rytov vs Maxwell co 0.659%는 scalar exact vs Maxwell co 0.648%와 거의 같다. 즉 잔차는 scalar 모델 차이이고 Rytov 절단은 사라졌다. Δn을 100배 줄이면 scalar 기준 Rytov 오차는 116배, Maxwell co 기준은 16배 감소한다. 크기 실험의 s=0.01과 이 t=0.01은 중심선 위상 δ≈0.021 rad로 같지만 Rytov 오차는 0.72% vs 0.09%, Maxwell co 기준 30% vs 0.66%로 전혀 다르다. **같은 δ라도 크기와 산란 각도 분포가 다르면 결과가 다르다.** 현재 형상과 넓은 NA를 유지하며 scalar Rytov를 쓰려면 Δn을 낮추는 쪽이 맞고, 이 t=0.01 조건이 inverse의 입력이 되었다.

추가로 t=0.01, θ=0–20°에서 scalar Rytov와 Maxwell exact의 xz magnitude·phase 5×4 그림을 저장했다(`2/low_contrast_xz_5x4.png`, phase color scale ±0.026 rad).

### 3.11 FFT·sampling 관련 아이디어와 검증

- **Nyquist guard**: Rytov 위상은 unpupilled propagating hemisphere 전체에서 만들어야 하므로 detector NA와 무관하게 dx < λ/(2n_m)=0.199 μm가 필요하다. 작은 detector NA만 보고 voxel pitch를 키울 수 없다.
- **Far amplitude → detector spectrum**: E_s ≈ A(k̂_s) e^{ikr}/r일 때 S(kx,ky;z) = 2πi A e^{i kz z}/kz. 이 변환과 Parseval로 exact 해와 FFT forward를 같은 detector spectrum에서 비교했다. r과 2r에서 추출한 A의 일치 ≤2.1e−6.
- **연속 형상 Fourier transform 대조군**: F(q) = f·V·3(sinχ−χcosχ)/χ³, χ²=b²(qx²+qy²)+a²qz². voxel 경계와 3D 보간 오차를 제거한 채 같은 1차 Born/Rytov 식을 평가한다. 이것이 없으면 voxel 오차와 근사 오차를 분리할 수 없다.
- **Padding/보간**: padding 2 → 3에서 voxel-연속 차이 4.5% → 2.0%(oblate). inverse의 matched operator에서는 padding 2/3/4에서 1.7% 대비 6.33% → 3.14% → 1.69%.
- **Window/sampling**: 25.6 → 51.2 → 102.4 μm window에서 spectrum 변화 0.074% → 0.034%, dx 0.1 → 0.05에서 0.0001%. 각도 coarsening(NA step 0.05→0.1, φ 45°→90°) 영향 ≤0.05 pp.
- **Evanescent 제외**: 모든 detector forward는 kz>0 propagating mode만 포함한다. z=5.2 μm 결과도 evanescent near-field 검증이 아니다.
- **Ray 위상 지표의 한계**: δ = k0·Δn·2a는 정렬 지표일 뿐 Rytov 정확도의 보편적 임계값이 아니다(3.3절 geometry control, 3.10절 크기/대비 비교).

---

## 4. 축 3: inverse — Rytov 데이터의 굴절률 복원 (`inverse/`)

### 4.1 v1: linear ODT (2026-09-12)

미지수 χ = n²/n_m² − 1, f = k_m²χ. 검출 전체장 U와 reference u0에서 Born은 U−u0, Rytov는 u0·log(U/u0)를 선형화 데이터로 쓰고, 2D detector FFT × dx²에 −2i kz e^{−i kz z_det}/k_m²를 곱해 χ의 Fourier sample을 얻는다(Fourier diffraction theorem). 각 sample은 Ewald 위치 q = k_s − k_i에 놓인다.

- `odt_prepare.m`: Ewald 위치를 **nearest-bin**으로 Cartesian Fourier grid에 놓고 중복은 평균. 실수 χ 가정으로 Hermitian 켤레 bin을 보충(`measured_mask`로 구분). Rytov 위상은 x/y 두 순서 unwrap 일치와 near-zero border로 2π branch를 정한다.
- `odt_reconstruct.m`: `direct`(zero-fill + IFFT), `gp`(Gerchberg–Papoulis, 측정 bin 대입 ↔ 비음수 투영), `tv`(비음수 isotropic 3D TV, Split Bregman; Goldstein–Osher penalized, Lim 2015 Eq. 15의 외부 재주입 없음).
- 검증: 중앙 voxel의 Fourier sample = χ0·dx³ (2e−12), displaced point scatterer (1e−11), 완전 데이터 역변환 (1e−12), 해석적 two-level TV 문제 (2e−6), 4.5 rad 위상 bump 복원.

데모(64³, dx=0.2, prolate 5×2.5, n=1.37, 49 angles, 잡음 없음), 굴절률 대비 상대 오차:

| 모델 | NA_det | NA_illum | direct | GP | TV |
|---|---:|---:|---:|---:|---:|
| Rytov | 0.7 | 0.2 | 0.607 | 0.602 | 0.634 |
| Born | 0.4 | 0.2 | 0.624 | 0.553 | 0.687 |
| Rytov | 0.7 | 0.5 | 0.464 | 0.781 | 0.520 |

TV가 수치적으로 수렴해도(745 iter) direct보다 나쁠 수 있음을 기록했다. 이후 forward의 Rytov pupil 순서 수정으로 이 v1의 dx=0.2 설정은 `RytovNyquistViolation`으로 멈추고, "pupil과 log가 교환한다"는 v1 test의 전제도 무효가 되었다.

### 4.2 v2: matched operator + exact Maxwell 데이터 (2026-09-13, `high_na_plan.md`)

조건: oblate a=3, b=5 μm, n_p=1.33567771866 (Δn/100), 81 방향(illumination NA≤0.5), detector NA=1.0, z=5 μm, detector 512×512 dx=0.1 μm, 복원 volume 128×128×80 dx=0.1 μm, Ewald sample 2,357,505개, 잡음 없음, 복원에 참 형상 미사용.

새 구성 요소:
- `odt_exact_data.m`: `spheroid-analytic-forward/2/case_3_reference.mat`의 validated Maxwell 해에서 co-polarized detector 전체장을 만든다. voxel truth나 inverse operator를 거치지 않은 **물리적 데이터**다. 검증 CSV(`case_3_reference_validation.csv`)와 run.log의 `PARAMETER_EXPERIMENT_PASS`를 확인한 뒤에만 로드한다.
- `odt_operator.m`: nearest-bin 대신 **padded unitary FFT + trilinear 보간과 그 정확한 transpose**. adjoint identity <1e−12, 기존 voxel forward와 <1e−11 일치. direct gridding은 padding과 무관하게 실제 복원 grid를 쓰고(초기 구현 2의 결함 수정), 밀도 정규화 + Hermitian 완성 + zero fill.
- `odt_reconstruct.m` v2: GP는 CGLS 최소 norm 데이터 보정 + 비음수, TV는 Lim Eq. 14 inner splitting + Eq. 15 residual 재주입, 정규방정식은 CG로.
- `run_inverse_experiment.m`: `source='exact'`(물리 데이터)와 `source='matched'`(참 χ를 operator에 통과시킨 선형 대조군)를 같은 solver로 돌려 **수치 오차와 모델 오차를 분리**한다. TV weight는 데이터 scale에서 유도한다.
- `check_inverse_sampling.m`: 연속 형상 vs 물리 데이터 vs voxel 예측 분리.

Operator 수렴 (`operator_convergence.csv`):

| padding | dx | 연속 형상 → 물리 데이터 | voxel → 연속 형상 | voxel → 물리 데이터 |
|---:|---:|---:|---:|---:|
| 2 | 0.1 | 0.570% | 6.33% | 6.23% |
| 3 | 0.1 | 0.570% | 3.14% | 3.06% |
| 4 | 0.1 | 0.570% | 1.69% | 1.67% |
| 3 | 0.08 | 0.570% | 3.05% | 2.97% |

0.570%는 1차 Rytov 모델 자체의 물리 오차(forward 실험의 0.659%와 일관)이고, 나머지는 trilinear 보간 오차다. detector 자체는 dx 0.05에서 2.9e−11, window 102.4 μm에서 1.4e−6 변화로 수렴했다.

**결과 스냅샷 (구현 3, padding 4, GP 100 iteration; 2026-09-13 01:43 기준)** — `exact_metrics.csv`, `linear_control_metrics.csv`:

| 항목 | exact, direct | exact, GP | matched, direct | matched, GP |
|---|---:|---:|---:|---:|
| 참 물체 → 물리 Rytov 데이터 불일치 | 1.675% | 1.675% | 1.675% | 1.675% |
| 굴절률 대비 상대 오차 ‖n−n_true‖/‖n_true−n_m‖ | 46.91% | 36.52% | 46.63% | 35.84% |
| 물체 내부 / 배경 오차 | 41.0% / 22.8% | 33.2% / 15.3% | 40.9% / 22.4% | 32.0% / 16.2% |
| 물체 평균 n bias (참값 1.33567772) | −10.4% | −5.9% | −11.5% | −8.3% |
| 축 방향 half-contrast span (참값 5.9 μm) | 4.3 μm | 8.0 μm | 4.3 μm | 8.0 μm |
| IoU | 0.868 | 0.864 | 0.868 | 0.874 |
| Ewald 데이터 잔차 | 18.08% | 1.19% | 18.04% | 0.15% |
| 물리 산란장 잔차 | 18.13% | 1.19% | 19.00% | 1.70% |
| 수렴 | direct | 미수렴 (100 iter, 270 s) | direct | 미수렴 (100 iter, 160 s) |

해석:
- exact와 matched의 결과가 거의 같으므로 direct의 47% 오차는 Rytov 모델 오차(1.7%)가 아니라 **선형 역문제의 missing cone**(illumination NA 0.5, 고정 +z detector)에서 온다. 축 방향 두께 과소평가(4.3 vs 5.9 μm)와 굴절률 과소평가(−10%)가 전형적인 missing-cone 증상이다.
- GP(CGLS 데이터 보정 + 비음수)는 데이터 잔차를 18% → 1.2%로 거의 없애지만 굴절률 오차는 36.5%에 머문다. 축 방향 span은 이제 8.0 μm로 과대평가된다. **데이터를 잘 맞추는 것과 굴절률을 맞추는 것은 missing-cone 문제에서 다른 일**이다.
- matched GP의 물리 산란장 잔차 1.70%는 Ewald 잔차 0.15%보다 크다. 이 차이가 곧 trilinear 보간 + Rytov 모델 오차(1.675%)다.

**진행 중인 실행**: `linear_control`의 TV(outer 1/5 완료, 데이터 잔차 0.0048), exact TV, `padding8` 대조 실행이 이 보고서 작성 시점에 MATLAB 프로세스 3개에서 돌고 있다. TV 결과, 그림, `inverse/results/report.md`는 아직 없다. 초기 구현 2(padding 2)의 provisional 결과(direct 45.9%, GP 52.0%, TV 38.6%, 모두 미수렴)는 gridding grid가 padding에 따라 바뀌는 결함이 있어 최종 결과로 쓰지 않는다(`results/provisional/README.md`).

**v2 검증 상태** (`verification.md`, 2026-09-13): inverse test suite 전체 통과(물리적 pupil/log 처리, FFT 부호·scale, 복소 adjoint, full-FFT 동치, direct grid의 padding 불변성, 해석적 TV endpoint, Lim 재주입, GP 독립 pseudoinverse 비교, stationarity false-stop 회귀, oblique detector energy, voxel forward 일치). Code Analyzer 12개 파일 issue 0. 독립 read-only 리뷰로 TV stationarity 오류 은폐, GP test 부족, cache/label 검증, exact-data provenance 검사를 수정했다. `README.md`는 v2 기준으로 갱신되었고, `high_na_plan.md` 체크리스트는 experiment/analysis 항목이 아직 미완이다.

---

## 5. 핵심 아이디어와 교훈 정리

1. **오차 분리의 도구가 결론을 바꾼다.** Maxwell exact만 있으면 25%가 scalar 탓인지 Rytov 탓인지 알 수 없다. scalar 경계값 exact(0.25% 차이), 연속 형상 transform(voxel 오차 제거), vector Born 대조군(작은 물체에서 scalar 모델 붕괴 판별), matched 선형 대조군(inverse의 모델/수치 분리)을 각각 만들었다.
2. **Rytov는 "큰 위상"이 아니라 "위상 기울기"에 깨진다.** 생략항 Q=∇ψ·∇ψ는 중앙과 가장자리의 위상 차이가 좁은 횡방향 거리에서 생길 때, 내부 반사파가 축에 모여 전진파와 간섭할 때, 그리고 field zero 근처에서 크다. 같은 4.09 rad이라도 평판 3%, 구 13%, prolate 24%.
3. **NA**: (a) Rytov의 detector pupil은 지수화 뒤에. (b) scalar 축약 정확도는 성분(co vs full)과 NA에 따라 다르고, 비정상입사에서는 작은 detector NA가 오히려 나쁘다. (c) illumination NA 0.5까지 Rytov 오차는 완만히 증가(9.9 → 11.2%). (d) 좁은 illumination NA는 inverse에서 missing cone으로 되돌아온다.
4. **FFT**: voxel/Ewald 보간 오차(padding 2에서 4–6%)는 근사 오차(10%)와 같은 자릿수이므로 반드시 수렴 검사가 필요하다. inverse에서는 nearest-bin gridding 대신 정확한 adjoint를 가진 trilinear operator를 써야 GP/TV의 데이터 항이 의미를 가진다.
5. **평가 지표를 고정하라.** 산란장 vs 전체장, 복소 vs 크기, pupil 유무, 평균 방식(가중 평균 vs pooled L2)에 따라 같은 계산이 25%, 10%, 8%, 5%로 보인다.
6. **크기 축소와 대비 축소는 다르다.** 같은 δ에서도 작은 물체는 scalar 모델을 깨고(30%), 낮은 Δn은 Rytov를 정확하게 만든다(0.09%). 넓은 NA를 유지한 채 scalar Rytov를 쓰려면 Δn을 낮추는 쪽이 맞다.
7. **2차 Rytov는 부분적 개선일 뿐** 먼 검출면에서 발산할 수 있다.
8. **Fail-closed 수치**: 검증 실패를 경고로 덮지 않고 error로 낸다. 이 정책 덕분에 좌표 roundoff(|η|>1), cutoff 부족((20,10) → (40,14)), NA 순서 버그 같은 실제 오류가 조용히 지나가지 않았다.

---

## 6. 저장소 상태와 재현

Git: `main`, HEAD `51094f1`. commit된 것은 spheroid solver(10 suite)와 scalar-polarization test뿐이다. `forward/`, `inverse/`, `spheroid-analytic-forward/1`, `2`, `ref/background`, `ref/tomocube-rytov`는 untracked이고 `sph_coords.m`(η clip), `test_sph_coords.m`, `spheroid-analytic-forward/history.md`는 modified다.

MATLAB R2024a를 `-nojvm -nodesktop -nosplash -r` + session-local `restoredefaultpath`로 실행했다. 수치 계산은 MATLAB에서 완료하고, 저장한 MAT를 Python(NumPy/SciPy/Matplotlib, `/tmp/mie-oblate-plot-env`)으로 그렸다.

| 목적 | 명령 (저장소 root, MATLAB) |
|---|---|
| exact solver 전체 test | `run('spheroid-analytic-forward/tests/run_all.m')` (약 2.5시간 이상) |
| forward 연속 형상 test | `run('forward/tests/run_all.m')` |
| voxel NA 순서 회귀 + prolate 15-case | `addpath('forward/tests'); test_born_rytov_voxel; analyze_rytov_spheroid;` |
| prolate 원인 진단 | `addpath('forward/exp','spheroid-analytic-forward'); test_scalar_prolate; diagnose_rytov(512); diagnose_rytov_term(.05); diagnose_rytov_geometry; locate_field_zero;` |
| oblate 81-direction sweep, xz, Q | `addpath('forward/exp'); oblate_na_sweep(stage)`, stage = check, reference, sweep, fieldcheck, report, xz, xz_validate 순서 |
| 크기 / 대비 실험 | `addpath('spheroid-analytic-forward/1'); run_size_experiment` / `addpath('spheroid-analytic-forward/2'); run_contrast_experiment` |
| inverse v2 | `run('inverse/tests/run_all.m'); addpath('inverse'); s=struct('padding',4,'gp_iter',100,'tv_inner',200,'cg_max_iter',100); run_inverse_experiment('exact',s); s.source='matched'; run_inverse_experiment('linear_control',s); check_inverse_sampling` |

주요 문서: `spheroid-analytic-forward/plan.md`·`history.md`, `forward/plan.md`·`history.md`, `forward/tests/rytov_analysis.md`, `forward/exp/report.md`·`oblate_report.md`·`oblate_xz_report.md`·`oblate_Q_mechanism.md`, `spheroid-analytic-forward/1/readme.md`·`2/readme.md`, `inverse/plan.md`·`high_na_plan.md`·`README.md`·`verification.md`.

---

## 7. 미완료 항목과 다음 단계

- **inverse v2의 TV 완주와 padding 8 대조**(진행 중), 그림과 한국어 결과 보고서(`inverse/results/report.md`). missing cone에 대해 GP는 데이터 잔차만 줄이고 굴절률 오차 36%와 축 방향 과대평가를 남겼다. TV(Lim 재주입)가 이를 얼마나 개선하는지가 이 프로젝트의 마지막 질문이다.
- TV weight·iteration 민감도, FOV/voxel pitch 민감도 실험(계획에는 있으나 미실행).
- `forward/tests/test_born_rytov.m:178`의 `InvalidTolerance` 기대와 현재 `born_rytov_spheroid.m` 검증의 불일치 해소.
- exact solver 11 suite 통합 실행과 537-case stress campaign.
- Q의 공간 source를 detector까지 전파해 위치별 오차 기여를 계산하는 분석, 반사 차수(Debye) 분해.
- 2차 Rytov의 production 구현 여부 판단(현재는 진단용).
- 실제 objective의 편광 전달·수차·analyzer, 유한 camera window, evanescent 성분, focused illumination.
- 균질 단일 spheroid를 넘어선 layered/heterogeneous 세포 모델용 reference solver.
- untracked 작업(`forward/`, `inverse/`, 실험 폴더)의 commit.
