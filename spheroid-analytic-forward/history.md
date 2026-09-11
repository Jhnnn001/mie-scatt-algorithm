# `spheroid-analytic-forward` 개발 및 검증 기록

작성 기준일: 2026-09-12

## 1. 이 문서의 목적

이 문서는 균질 spheroid의 전자기 산란장을 계산하는 MATLAB 코드를 어떤 근거로 설계했고, 구현 중 어떤 실패를 겪었으며, 무엇을 고쳐 최종 검증에 이르렀는지를 기록한다. 현재 구현은 spheroidal-coordinate separation of variables와 표면 Galerkin matching으로 Maxwell 방정식의 full-vector 해를 수치적으로 구한다. 여기서 "exact"는 symbolic closed-form 해가 아니라, truncation·quadrature·ODE·경계조건 검사를 통과한 수치적 full-wave reference를 뜻한다.

이 프로젝트의 두 목적은 다음과 같다.

1. prolate·oblate spheroid에 대한 복소 전기장 `(Ex,Ey,Ez)`를 임의의 입사 방향과 편광, 내부·외부 관측점에서 계산한다.
2. 그 full-vector 해를 기준으로 Born, Rytov, multi-slice BPM 같은 근사법의 유효 범위를 비교한다.

두 번째 목적은 일부 완료되었다. fixed-polarization test에 이어 biological prolate baseline에서 Born과 first Rytov의 복소 scattered-field error를 full-vector `Ex`와 비교했고, 제한된 experimental scalar prolate boundary-value reference로 vector-reduction error와 Rytov-truncation error도 분리했다. 그러나 geometry와 물리 parameter 전반의 systematic sweep가 없으므로 일반적인 유효 범위는 아직 정해지지 않았다.

상세 수학 명세와 현재 인터페이스는 [plan.md](plan.md)에 있다. 이 문서는 명세를 반복하는 대신, 명세가 왜 현재 형태가 되었는지와 어떤 증거가 남아 있는지를 설명한다.

## 2. 문헌, 물리 convention, 범위

주요 문헌은 다음과 같다.

- Asano and Yamamoto, *Light Scattering by a Spheroidal Particle*, Applied Optics 14, 29 (1975): [로컬 PDF](../ref/mie-scatt/ao-14-1-29.pdf)
- Barton, *Internal, near-surface, and scattered electromagnetic fields for a layered spheroid with arbitrary illumination*, Applied Optics 40, 3598 (2001): [로컬 PDF](../ref/mie-scatt/ao-40-21-3598.pdf)
- DLMF Chapter 30과 Flammer (1957): spheroidal eigenvalue, angular function, radial function의 convention과 recurrence 확인
- Bohren and Huffman 및 Wiscombe Mie benchmark: sphere fallback과 광학정리 convention 확인

모든 구현은 시간 convention `exp(-i omega t)`를 사용한다. 진공 파장은 `lambda`, 배경과 입자의 굴절률은 각각 `n_m`, `n_p`이며, `n_p`에는 흡수를 나타내는 양의 허수부를 허용한다. spheroid의 대칭축은 전역 z축이고, `a`는 z 방향 semi-axis, `b`는 xy 평면의 semi-axis이다. 모든 길이는 같은 단위를 사용하며 내부 특수함수 계산에서는 semifocal length로 무차원화한다.

## 3. 저장소 이력

| 시각 (KST) | commit | 의미 |
|---|---|---|
| 2026-09-11 23:22 | `97618ee4544c2f56889682a3acd2eb7e768defdc` | `spheroid-analytic-forward/` 디렉터리 생성 |
| 2026-09-12 04:36 | `5cecb3f8ef326150d77155a93071f9a7a4d1e483` | 초기 plan, 안정화 좌표, 평면파, Gauss–Legendre 구현 |
| 2026-09-12 14:09 | `569307bbee7f67b14ee07eb97c649a96681e3dd4` | full spheroidal solver와 당시 10개 검증 suite 구현 및 push |
| 2026-09-12 15:38–17:49 | uncommitted | biological scalar-polarization 6-case characterization |
| 2026-09-12 18:45–18:56 | uncommitted | detector-NA Rytov operator 수정과 15-case exact-vector comparison |
| 2026-09-12 19:10–19:18 | uncommitted | experimental scalar prolate boundary-value reference와 초기 Rytov truncation 진단 |
| 2026-09-12 19:18–19:24 | uncommitted | grid refinement, sphere geometry control, field-zero localization, final experiment rerun |

`569307b`는 작성 시점에 `main`과 `origin/main`이 공유하는 solver commit이다. 그 뒤 추가된 biological scalar-polarization test와 이에 맞춘 [plan.md](plan.md), [tests/run_all.m](tests/run_all.m) 변경은 아직 commit되지 않았다. `forward/`의 Born/Rytov 구현·분석도 untracked 상태다. 따라서 "commit된 10개 suite", "별도로 완료된 scalar-polarization suite", "별도로 완료된 Born/Rytov 분석"을 하나의 통합 test run이나 commit된 결과로 합쳐 말하면 안 된다.

## 4. 초기 설계 검토에서 바뀐 핵심 내용

초기 plan은 수학적 방향은 맞았지만, 극좌표 특이점과 고차 특수함수의 부동소수점 거동을 과소평가했다. 세 차례의 검토를 거쳐 다음 항목을 명세에 반영했다.

### 4.1 좌표 변환

초기 oblate 변환은 `z=0`에서 `sign(z)=0`을 그대로 사용하여 focal disk의 모든 점을 `eta=0`으로 잘못 보낼 수 있었다. `sign(0):=+1`로 sheet를 고정했고, radial·angular parity 때문에 반대 sheet를 선택해도 같은 물리장이 나와야 한다는 점을 확인했다.

더 큰 문제는 prolate 축과 초점선 근처에서 `xi^2-1`, `1-eta^2`를 큰 수의 차로 계산하면 유효숫자를 잃는다는 것이었다. 해결책은 `alpha=xi^2-1` 또는 `xi^2+1`, `beta=1-eta^2`, `D`를 상쇄 없는 식으로 직접 계산하여 [sph_coords.m](sph_coords.m)이 반환하고, 모든 후속 함수가 `xi`, `eta`에서 이 작은 양을 재계산하지 않도록 한 것이다.

초점과 focal ring은 좌표계 자체의 singular set이므로 제거할 수 없다. 그 주변 입력은 무차원 거리 `delta=1e-8` 밖으로 옮기되, 방위각과 접근 side를 보존하고 이동 후에도 singular set에서 최소 `delta`만큼 떨어지도록 구현했다. 축, prolate focal segment, oblate disk 자체는 regular limit를 사용하므로 임의로 이동시키지 않는다.

### 4.2 outgoing radial function

초기 설계는 spherical Neumann series로 `R^(2)`를 계산하고 `R^(3)=R^(1)+iR^(2)`를 구성하려 했다. 그러나 항 비율이 점근적으로 `1/xi^2`이므로 `xi`가 1에 가깝거나 oblate에서 `xi<1`이면 수렴이 지나치게 느리거나 발산하며, 그 전에 `y_n` overflow와 계수 underflow가 발생했다. 작은 expansion coefficient를 기준으로 항을 버리는 방식도 Neumann series에는 안전하지 않았다.

중간안은 큰 `xi`의 series와 inward ODE를 중첩 구간에서 연결하는 것이었다. 실제 구현에서는 canonical normalization 자체가 boundary matching에 필요하지 않다는 점을 이용해 더 단순하고 안정적인 방향을 택했다. 최종 exterior basis는 mode마다 `R_out(xi0)=1`로 임의 정규화하고, outgoing inverse-`xi` asymptotic seed에서 Riccati/log-amplitude ODE를 안쪽으로 적분한다. Galerkin coefficient가 이 임의 scale을 흡수하므로 물리장은 변하지 않는다.

이 선택 때문에 canonical Wronskian 값은 더 이상 올바른 검사가 아니다. 대신 서로 다른 두 시작점에서 적분한 normalized log derivative와 field의 일치, 그리고 ODE·asymptotic backward error를 독립적으로 검사한다. 검사가 실패하면 경고로 덮지 않고 해당 basis를 거부한다.

### 4.3 수렴과 검증 정책

표면에서 Galerkin system을 풀었다는 사실만으로 수렴을 주장할 수 없다. 최종 정책은 incident azimuthal spectrum tail, `Q_phi` 독립성, 상위 세 `l` shell의 실제 field contribution, 독립 surface grid의 boundary residual, `Q_theta` 재계산, `(L+5,M+2)` 재계산, 특수함수 진단을 모두 통과해야 `validated=true`를 허용한다.

표면 적분은 `eta`에서 직접 Gauss–Legendre를 적용하지 않고 `eta=cos(theta)`로 바꿨다. 홀수 azimuthal order에서는 `sqrt(1-eta^2)` endpoint behavior 때문에 `eta` 적분의 수렴이 대수적으로 느려질 수 있지만, `theta` 적분에서는 integrand가 매끄럽다. `sin(theta)` Jacobian을 명시적으로 포함한다.

복소 `n_p`에서는 eigenvalue label이 바뀔 수 있다. 완전한 retained subspace 안의 순열은 Galerkin 해를 바꾸지 않으므로 모든 mode를 `c=0`부터 추적하지 않았다. 대신 eigenvector matrix condition, coefficient가 아닌 실제 modal-field tail, 증가한 cutoff에서 재구성 field와 observable의 변화를 검사한다.

### 4.4 angular derivative

초기 reduced Legendre derivative 공식은 비정규화 Ferrers 함수의 계수를 정규직교 기저에 그대로 적용한 오류가 있었다. 정규직교 reduced basis에서는 order를 하나 또는 둘 올릴 때 degree-dependent square-root factor가 필요하다. 이를 term별로 적용하고, raw `S'`, `S''` 대신 vector field가 실제로 요구하는 regular weighted combinations `W0`, `W1`, `W2`, `Wm`을 직접 계산하도록 바꿨다.

특히 `m=1`에서는 `S''`뿐 아니라 `(1-eta^2)S''`도 endpoint에서 발산할 수 있다. 따라서 `(1-eta^2)^(3/2) S''`까지 metric factor와 결합한 뒤 평가해야 axis에서 유한한 vector component를 얻는다.

### 4.5 지원 범위와 물리 test

단순히 `c<=200`이라고 쓰는 대신 외부와 내부 size parameter를 분리했다. 현재 vector-solver runtime-validation envelope의 주요 parameter bound는 `x_ext=k_m max(a,b)`와 `x_int=|k_p|max(a,b)`가 각각 250 이하, aspect ratio `a/b`가 1/4–4, `Re(n_p/n_m)`가 0.5–2, `Im(n_p)/Re(n_p)`가 0–0.1인 범위다. 이 범위 안이라는 사실만으로 검증되는 것은 아니며, 매 실행의 수치 검사를 모두 통과해야 한다.

zero contrast, Rayleigh limit, lossless optical theorem, lossy Poynting/volume absorption balance, 내부·외부 boundary trace, Mie sphere와 near-sphere limit, Barton benchmark를 필수 test로 추가했다.

## 5. 구현의 실제 진행

### 5.1 04:12–04:39 — 기하, quadrature, 입사장

[gauss_legendre.m](gauss_legendre.m), [sph_coords.m](sph_coords.m), [incident_plane_wave.m](incident_plane_wave.m)을 먼저 구현했다. 이 단계에서 arbitrary incidence direction과 complex polarization을 정의하고, 입사 전기장과 자기장이 `H=n_m khat cross E`를 만족하는지 검사했다.

좌표 test는 처음부터 여러 실제 오류를 드러냈다. singular-neighborhood mask가 너무 넓어 정상점을 옮겼고, subnormal 좌표의 제곱이 underflow했으며, focal ring에 접선 방향으로 접근할 때 상쇄가 발생했고, 첫 nudge가 반올림 때문에 정확히 ring에 남는 경우도 있었다. 초기 단계에서는 strict mask, subnormal-safe reconstruction, tangential-ring 처리와 두 번째 nudge로 이를 막았다. 방위각·side를 보존하는 guaranteed-distance projection과 overflow-safe semifocal length [sph_semifocal.m](sph_semifocal.m)은 뒤의 hardening 단계에서 추가되었다.

### 5.2 04:39–04:57 — sphere/Mie fallback

구가 spheroidal coordinate에서 semifocal length 0으로 퇴화하므로, 구에는 별도의 [mie_field.m](mie_field.m)과 [mie_efficiencies.m](mie_efficiencies.m)을 사용했다. near-equal axes는 floating-point spacing 기준으로만 sphere route에 보내며, 임의의 넓은 near-sphere 근사는 사용하지 않는다.

초기 boundary test와 zero-contrast test는 vector spherical harmonic의 radial `N` 성분에 `pi_n` factor가 빠진 오류를 발견했다. 이 오류는 tangential field 일부가 그럴듯해 보여도 normal displacement continuity와 내부 incident-field reconstruction을 깨뜨렸다. factor를 복구한 뒤 tiny-size guard, internal-size truncation, outgoing phase sign, interior-axis test를 추가했다.

### 5.3 04:52–05:18 — eigenproblem, angular function, Bessel, regular radial function

[sph_eigen.m](sph_eigen.m)은 even·odd parity별 complex-symmetric tridiagonal eigenproblem을 푼다. angular expansion은 Condon–Shortley phase가 없는 orthonormal Ferrers basis를 사용하며, 복소 `c`에서도 bilinear normalization을 유지한다. truncation을 1.5배 늘려 retained eigenvector가 안정적인지 재검사한다.

[sph_angular.m](sph_angular.m)은 endpoint singular division을 피하기 위해 reduced Ferrers polynomial과 weighted derivative를 계산한다. [sph_bessel.m](sph_bessel.m)은 spherical Bessel 함수와 도함수를 order-by-argument 배열로 제공하고, [sph_bessel_T.m](sph_bessel_T.m)은 oblate `xi<1`에서 `T_r(x)=j_(mu+r)(x)/x^mu`를 log-safe 또는 small-`x` series로 계산한다.

여기서 MATLAB의 `besselj`가 order vector와 argument vector를 원하는 outer-product 형태로 자동 broadcasting하지 않는다는 문제가 나타났다. 또한 고차에서 개별 표현은 `0*Inf`가 되어 Wronskian을 계산할 수 없고, ordinary `j_n`이 underflow해도 scaled `T_r`는 유한할 수 있으며, raw expansion coefficient tail이 작다고 실제 radial contribution이 작다는 보장도 없었다. 배열을 명시적으로 구성하고, scaled representation과 실제 term contribution 기준을 사용하도록 수정했다.

prolate `c=50`의 일부 regular radial mode는 normalization sum condition number가 `1.76e15–7.08e15`에 달했다. 식의 phase나 factorial을 임의로 바꾸는 것으로 해결할 문제가 아니었고, parity/Frobenius initial-value evaluation을 fallback으로 도입했다. 반대로 high-order complex oblate 예 `c=200+20i, L=234`는 표현 민감성이 남아 있으므로 조용히 부정확한 값을 내지 않고 `sph_eigen:EigenSensitivity` error를 낸다.

### 5.4 05:19–06:13 — vector spheroidal wave와 outgoing ODE

[sph_vecwave.m](sph_vecwave.m)은 scalar separated solution에서 regular vector wave `M`과 `N`을 구성한다. 모든 식을 `W0–W2`, `V0–V2`, `alpha`, `beta`, `D`로 다시 써 axis와 focal disk에서 불필요한 `0/0`을 제거했다. Cartesian finite difference로 `curl M=cN`, `curl N=cM`, `div M=div N=0`을 prolate·oblate, 양·음 `m`, axis·disk·generic point에서 확인했다.

이 단계에서 canonical `R^(2)` series와 splice 접근을 폐기하고 [sph_radial_ode.m](sph_radial_ode.m), [sph_radial_ode_eval.m](sph_radial_ode_eval.m)의 임의 정규화 outgoing basis로 전환했다. 상태 변수는 큰 oscillatory term을 제거한 `w=R'/R-ic+1/xi`와 phase-removed log amplitude를 사용한다. 80항 이하의 inverse-`xi` recurrence로 seed를 만들고, 두 독립 시작점의 agreement가 `1e-10`보다 작아질 때까지 시작점을 바깥으로 옮긴다.

`c<1`에서는 asymptotic coefficient 자체가 overflow하지 않도록 `b_n=a_n c^n`을 저장한다. `ode113` dense output을 surface부터 시작점까지 사용하고, 그 밖에서는 같은 asymptotic expansion으로 이어간다. oblate `xi0=0`은 log normalization이 singular하므로 의도적으로 지원하지 않으며, 현재 검증 범위의 oblate 하한은 그보다 큰 `xi0`를 가진다.

### 5.5 05:36–07:34 — Galerkin solver와 public API

[spheroid_solve.m](spheroid_solve.m)은 표면에서 tangential `E`, tangential `H`를 incident field와 일치시키는 Galerkin system을 각 signed azimuthal order `m`별로 조립한다. 내부 regular coefficients와 외부 outgoing coefficients를 동시에 풀고, normal `n^2 E` continuity는 독립 검증으로 남긴다.

초기 matrix는 `rcond`가 `1e-20` 부근까지 떨어졌다. arbitrary-normalized basis의 column-scale disparity는 column equilibration으로 처리했지만 낮은 `rcond`가 남았고, Rayleigh 분석은 일부 활성 system의 구조적 rank deficiency도 확인했다. 작은 RHS만 보고 해당 mode를 0으로 만드는 시도는 Rayleigh case를 망가뜨렸으므로 폐기했다. 최종 구현은 axial incidence에서 해석적으로 존재하지 않는 sector만 정확히 0으로 두고, 활성 matrix가 rank-deficient이면 `lsqminnorm`을 사용하며 rank와 solve method를 기록한다.

[spheroid_eval.m](spheroid_eval.m)은 저장된 solution을 내부·외부 점에서 평가한다. `side='auto'`는 Cartesian ellipsoid equation으로 boundary를 판정하고, `side='in'|'out'`은 boundary trace 검사용으로 제한한다. `'scattered'`는 외부에서 작은 두 수의 subtraction을 하지 않고 modal expansion을 직접 반환한다. malformed active sector, 손상된 ODE payload, nonfinite coefficient는 fail-closed error가 된다.

[spheroid_field.m](spheroid_field.m)은 입력 검증, sphere/prolate/oblate routing, solve, evaluation을 묶는 public entry point다. zero contrast, incident-only 요청, sphere에는 불필요한 spheroidal solve를 피하는 analytic fast route가 있다.

### 5.6 07:15–13:53 — post-implementation hardening

완성 직전 audit에서 oblate scaled derivative sign, regular-series vectorization, surface normal scaling, stored-solution validation, axial metadata, Mie tolerance metadata를 다시 고쳤다. 이때 [sph_semifocal.m](sph_semifocal.m)의 overflow-safe scale 계산과 방위각·side를 보존하는 guaranteed-distance singular projection도 최종 형태로 들어갔다. 길이를 공통으로 `1e-200` 또는 `1e200` 배 해도 동일한 무차원 물리가 나오는지, polarization amplitude를 같은 극단값으로 바꿔도 선형성이 유지되는지 검사했다.

stress report는 중간 종료 후 재개할 수 있도록 row count와 case identity를 검증하고 MAT 파일을 atomic하게 교체한 뒤 CSV를 파생하도록 작성했다. 하지만 537개 전체 case를 실제 solver로 실행한 것은 아니다. stub solver로 case definition과 report/resume mechanism만 smoke-test했다.

중간 full run은 Rayleigh outgoing verification과 Barton case 하나에서 실패했다. outgoing diagnostic의 roundoff floor, radial tolerance scaling, far-field sampled-peak 평가를 수정한 뒤 당시 10개 suite가 모두 통과했다. 이 상태가 commit `569307b`다.

## 6. 최종 solver 구조

| 계층 | 파일 | 역할 |
|---|---|---|
| Public API | [spheroid_field.m](spheroid_field.m) | 입력 검증, route 선택, solve와 field evaluation |
| Solution | [spheroid_solve.m](spheroid_solve.m) | Galerkin assembly, adaptive cutoff·quadrature, fail-closed validation |
| Evaluation | [spheroid_eval.m](spheroid_eval.m) | 내부·외부 total/scattered/incident `E,H` 계산 |
| Geometry | [sph_semifocal.m](sph_semifocal.m), [sph_coords.m](sph_coords.m) | overflow-safe scale과 cancellation-free spheroidal coordinates |
| Angular | [sph_eigen.m](sph_eigen.m), [sph_angular.m](sph_angular.m) | eigenvalue·Legendre coefficients와 regular weighted angular combinations |
| Radial | [sph_bessel.m](sph_bessel.m), [sph_bessel_T.m](sph_bessel_T.m), [sph_radial.m](sph_radial.m) | regular radial basis와 scaled small-`xi` path |
| Outgoing radial | [sph_radial_ode.m](sph_radial_ode.m), [sph_radial_ode_eval.m](sph_radial_ode_eval.m) | inverse-`xi` seed, inward Riccati/log ODE, far continuation |
| Vector basis | [sph_vecwave.m](sph_vecwave.m) | full-vector `M,N` spheroidal wave functions |
| Sphere route | [mie_field.m](mie_field.m), [mie_efficiencies.m](mie_efficiencies.m) | exact sphere fallback와 benchmark quantities |
| Illumination | [incident_plane_wave.m](incident_plane_wave.m) | arbitrary direction·complex polarization의 plane wave |
| Demonstration | [example_xz_slice.m](example_xz_slice.m) | Barton형 prolate case의 x–z field slice |

외부 radial basis의 절대 normalization은 물리량이 아니다. surface basis와 coefficient의 곱만 물리적이며, 이 관찰이 불안정한 canonical Neumann normalization을 제거하게 한 핵심 설계 결정이다.

## 7. solver 검증 결과

commit `569307b`에 대응하는 보존 로그는 `../../tmp/spheroid-analytic-forward/green/run_all_final10.log`이다. 이 상대경로는 저장소의 형제 임시 디렉터리 `/Users/jhlee/Documents/jhlee_forge/projects/tmp/`를 가리킨다.

당시 10개 suite의 결과는 다음과 같다.

- Gauss–Legendre exactness, plane-wave convention, spheroidal coordinates, stored-solution validation, unit/amplitude scaling: 통과
- angular·radial special functions와 vector identities: 통과
- Wiscombe Mie benchmark, Rayleigh limit, sphere boundary, near-sphere comparison: 통과
- lossless prolate energy mismatch: `9.89e-8`
- lossless oblate energy mismatch: `4.43e-8`
- lossy prolate와 oblate에서 Poynting surface flux와 volume-loss integral로 각각 구한 absorption이 `C_ext-C_sca`와 허용오차 `1e-5` 안에서 일치
- base와 refined volume quadrature 사이의 상대 변화: prolate `5.23e-16`, oblate `3.18e-15`
- 최대 boundary component residual: prolate normal `n^2 E`의 `1.86e-7`
- Barton far-field sampled maximum 상대오차: `7.92e-8`, `5.66e-7`, `3.02e-7`, `4.25e-7`; 기준 `1e-4`보다 작음
- MATLAB Code Analyzer: 30개 파일에서 issue 0개
- [example_xz_slice.m](example_xz_slice.m): 90,601개 finite field magnitude 계산, 범위 `0.00432103841629–5.97273629666`; 수치 범위의 보존 로그는 `../../tmp/spheroid-analytic-forward/green/example_xz_slice_final3.log`, 최종 재실행 성공 로그는 `../../tmp/spheroid-analytic-forward/green/example_xz_slice_final5.log`

두 느린 suite의 당시 실행시간은 `test_physics`가 7,520.810 s, `test_barton2001`이 1,594.542 s였다. 따라서 작은 code change마다 전체 suite를 반복하기는 비싸며, focused test를 먼저 실행하고 release 또는 기준점 갱신 때 전체 suite를 실행하는 편이 맞다.

`stress_envelope.m`은 runtime-validity domain의 537개 corner·interior·resonance case를 정의한다. 전체 numerical campaign을 완료한 결과 파일은 없으며, report format과 resume logic의 smoke 결과만 `../../tmp/spheroid-analytic-forward/green/stress_envelope_smoke_final.log`에 남아 있다. 그러므로 현재 envelope는 "입력 자격조건과 실행별 검증 계약"이지, 모든 점을 이미 경험적으로 인증했다는 뜻이 아니다.

## 8. scalar approximation 조사

### 8.1 처음 질문과 metric 수정

처음에는 x-polarized plane wave를 넣었을 때 `Ex`가, y와 z에서도 각각 같은 Cartesian component가 지배적인지를 확인하려 했다. 이 질문은 full-vector field를 한 scalar component로 축약할 수 있는지 확인하는 데 필요한 조건이지만 충분조건은 아니다.

초기 pilot은 `x_ext=1`, `n_p/n_m=1.001`의 매우 약한 contrast에서 여섯 shape/polarization case를 계산했다. 좁은 patch의 leakage는 약 `0.052110–0.052119`, 넓은 patch는 `0.408500–0.413274`였고 약 49.3 s가 걸렸다. 이 값과 시간은 작업 session의 console 관찰값이며 별도 보존 로그는 남아 있지 않다. 이 조건은 실제 cell scale과 거리가 멀었고, Cartesian component 크기만 비교하면 amplitude·phase 오차와 sampling geometry를 숨길 수 있었다.

검토 후 [tests/test_scalar_approximation.m](tests/test_scalar_approximation.m)은 unit incident polarization `e0`에 대한 co-polarized projection

`E_co = (E_sca dot conj(e0)) e0`

와 aggregate leakage

`epsilon_perp = norm(E_sca-E_co,'fro') / norm(E_sca,'fro')`

를 사용하도록 바뀌었다. `retained=norm(E_co,'fro')/norm(E_sca,'fro')`도 함께 기록하지만, 이는 샘플된 복소 field 배열의 norm 비율이다. solid-angle 적분 에너지 비율이나 모든 방향에서의 pointwise bound가 아니다.

### 8.2 biological baseline

최종 test의 idealized cell parameter는 다음과 같다.

| 항목 | 값 |
|---|---|
| Vacuum wavelength | `0.532 um` |
| Background | 20 °C distilled water, `n_m=1.335381534` |
| Cell | homogeneous, isotropic, nonabsorbing, `n_p=1.37` |
| Prolate semi-axes | `a=5 um`, `b=2.5 um` |
| Oblate semi-axes | `a=2.5 um`, `b=5 um` |
| External size parameter | `x_ext=78.858` |
| Internal size parameter | `x_int=80.902` |
| Solver tolerance | `1e-6`, failure policy `error` |

x와 y polarization은 `+z`로 진행하고, z polarization은 Maxwell transversality `E0 dot khat=0` 때문에 `+x`로 진행한다. 따라서 z 결과와 x/y 결과의 차이는 polarization 효과만이 아니라 spheroid symmetry axis에 대한 incidence orientation 변화도 포함한다. x와 y는 rotational symmetry 때문에 중복 검증에 가깝다.

관측 거리는 `1e4 max(a,b)`이고, 방향은 `khat+u e0+v(khat cross e0)`를 정규화해 만든다. 좁은 정사각 patch는 `|u|,|v|<=0.08`로 축 방향 최대 4.574°, corner 6.455°이며, 넓은 patch는 `|u|,|v|<=1`로 축 방향 45°, corner 54.736°다. 이는 gnomonic square sampling이며 원형 numerical-aperture cone이나 solid-angle-weighted quadrature가 아니다.

### 8.3 계산 결과

| Shape | Polarization | Narrow leakage | Narrow retained | Wide leakage | Wide retained |
|---|---:|---:|---:|---:|---:|
| Prolate | x | 0.033456 | 0.999440 | 0.145096 | 0.989418 |
| Prolate | y | 0.033456 | 0.999440 | 0.145096 | 0.989418 |
| Prolate | z | 0.019954 | 0.999801 | 0.126329 | 0.991988 |
| Oblate | x | 0.020852 | 0.999783 | 0.132063 | 0.991241 |
| Oblate | y | 0.020852 | 0.999783 | 0.132063 | 0.991241 |
| Oblate | z | 0.034708 | 0.999397 | 0.203252 | 0.979127 |

모든 narrow case는 운영상 threshold `epsilon_perp<0.1`을 통과했고, 모든 wide case는 `epsilon_perp>0.1`이었다. 따라서 이 한 parameter set의 좁은 forward patch에서는 산란장을 고정 입사 편광 방향으로 투영하는 축약이 aggregate norm 기준으로 작동하지만, 검사한 wide forward gnomonic patch에서는 같은 기준을 통과하지 못했다.

이 각도 의존성에는 산란체뿐 아니라 far-field transversality의 기하학도 들어 있다. 원거리 산란장은 각 관측 방향 `rhat`에 대해 `E_sca dot rhat=0`이어야 하므로, 관측 방향이 크게 변하면 하나의 고정된 lab-frame `e0`와 계속 평행할 수 없다. 따라서 narrow leakage가 작은 것은 forward geometry에서 자연스러운 결과이기도 하고, wide leakage가 큰 것만으로 scalar Helmholtz·Born·Rytov amplitude model이 실패했다고 결론낼 수도 없다.

예를 들어 prolate x/y의 narrow leakage `0.033456`은 discarded vector의 amplitude-norm이 전체의 3.35%라는 뜻이다. 직교 projection의 norm-squared 비율로 환산하면 약 0.112%지만, 이를 검출기 power error로 그대로 해석하려면 solid-angle weight, optical transfer, detector model이 추가로 필요하다.

여섯 case는 약 2시간 11분 걸렸다. x/y axial incidence에서는 symmetry로 `|m|=1` sector만 필요하지만, z polarization을 만들기 위한 `+x` side incidence는 모든 azimuthal sector를 활성화한다. 작업 중 solver metadata에서는 마지막 oblate z case가 대략 `M=108`, `L=109`까지 진행하여 최대 217개의 signed `m` system과 validation re-solve를 요구한 것으로 관찰되었다. 최종 scalar log에는 이 cutoff가 출력되지 않았으므로 해당 두 수는 로그로 재검증할 수 없지만, 모든 `m`이 활성화되는 구조가 마지막 계산이 유난히 느렸던 직접 원인이다.

완료 로그는 `/Users/jhlee/Documents/jhlee_forge/projects/tmp/spheroid-analytic-forward/biological_scalar.log`에 있다. 이 실행은 `test_scalar_approximation: PASS`로 끝났다. 이후 해당 test와 runner에 대한 정적 검사도 작업 session에서는 issue 없이 끝났지만, 그 출력은 별도 로그로 보존되지 않았다.

### 8.4 Born/Rytov complex-field comparison

polarization test 뒤에는 실제 scalar forward model의 amplitude와 phase를 포함하는 비교를 수행했다. 결과와 방법은 [forward/tests/rytov_analysis.md](../forward/tests/rytov_analysis.md), 원자료는 [rytov_analysis.csv](../forward/tests/rytov_analysis.csv), 실행 로그는 [rytov_analysis.log](../forward/tests/rytov_analysis.log)에 보존되어 있다. 이 분석은 uncommitted `forward/` tree에 있으므로 commit `569307b`의 일부는 아니다.

분석 과정에서 [born_rytov_voxel.m](../forward/born_rytov_voxel.m)의 detector-NA 적용 순서에 실제 오류가 발견되었다. 이전 구현은 linear Born field를 detector pupil로 먼저 제한한 뒤 그 field로 Rytov phase를 만들었고, exponentiation 뒤 생긴 out-of-pupil frequency를 최종 출력에 남겼다. 그러나 first Rytov phase는 objective 이전의 unpupilled propagating field에서 구성해야 하며, exponentiation 뒤 detector pupil을 마지막에 적용해야 한다.

수정된 순서는 다음과 같다.

1. object FFT를 forward-propagating Ewald hemisphere에서 sampling한다.
2. pupil을 적용하지 않은 linear field를 inverse transform하고 incident field로 나누어 first Rytov complex phase를 만든다.
3. `u0.*expm1(phase)`로 scattered Rytov field를 구성한다.
4. 다시 Fourier transform한 뒤 detector NA pupil을 적용하고 detector plane으로 inverse transform한다.

Born field에는 최종 linear pupil만 적용한다. 두 출력은 같은 detector pupil로 band-limit된다. 정확한 순서를 고정한 회귀 test [test_born_rytov_voxel.m](../forward/tests/test_born_rytov_voxel.m)은 detector NA를 바꿨을 때 작은 pupil 안의 Born과 Rytov spectrum 변화가 각각 `2.194e-16`, `2.112e-16`, Rytov의 out-of-pupil leakage가 `1.856e-16`임을 확인했다. zero contrast, tilted illumination, independently summed translated point source, Nyquist와 object-spectrum bandwidth rejection도 통과했다.

baseline 비교는 prolate `a=5 um`, `b=2.5 um`, `lambda=0.532 um`, `n_m=1.335381534`, `n_p=1.37`, `+z` incidence, lab x polarization을 사용했다. detector는 forward tip보다 2 um 뒤인 `z=7 um`, `NA=0.10`, voxel grid는 `128^3`, `dx=0.1 um`, padding factor 2, field of view는 `12.8 um`였다. centerline ray phase-delay estimate는 `4.08861 rad`다.

reference는 [spheroid_solve.m](spheroid_solve.m)와 [spheroid_eval.m](spheroid_eval.m)의 validated scattered field를 `Ex`에 projection한 값이다. 각 detector spatial frequency의 outgoing direction에서 반경 `r`과 `2r`, `r=1e6 max(a,b)`로 far amplitude를 추출했고, 두 추출값은 15개 analysis row 전체에서 최대 `2.12e-6` 상대차로 일치했다. 사용된 exact spheroidal solution은 모두 `tol=1e-6` runtime validation을 통과했다.

reported error는 모두 detector pupil 안의 complex scattered-field relative L2 error

`norm(model-exact)/norm(exact)`

다. intensity, total field, image reconstruction error가 아니다. baseline 결과는 다음과 같다.

| Model | Relative complex scattered-field error |
|---|---:|
| Born | 197.56% |
| 수정 전 returned Rytov | 74.40% |
| 수정 전 Rytov에 final pupil만 추가 | 26.21% |
| 수정된 voxel Rytov | 25.74% |
| analytic continuous-shape first Rytov | 24.80% |

74.40%에서 25.74%로 보이는 큰 변화의 대부분은 out-of-pupil frequency를 제거한 효과다. premature pupil 자체를 고친 순수 in-pupil 변화는 26.21%에서 25.74%다. field of view가 25.6 um인 control에서는 final-pupil-only error가 41.98%, 올바른 순서가 25.33%여서 window에 따라 premature cutoff의 영향이 커질 수 있다.

continuous-shape control은 homogeneous spheroid의 analytic Fourier transform을 사용하여 voxel boundary와 3-D interpolation을 제거한다. 이는 exact scalar scattering solution이 아니라 동일한 first Rytov approximation의 독립적인 numerical control이다. baseline에서 voxel과 continuous Rytov field 차이는 padding 2에서 3.41%, padding 4에서 1.44%였고, continuous result는 field of view 51.2–102.4 um에서 24.30%–24.34%에 머물렀다. 따라서 약 24%의 residual을 voxelization이나 finite window만으로 설명할 수 없다.

물리 control도 같은 결론을 보강했다. relative index contrast를 `0.001`로 낮추면 corrected voxel Rytov error는 1.21%, fine continuous result는 1.04%로 줄었다. 반면 baseline contrast에서 wavelength를 633 nm로 바꾼 fine continuous error는 22.84%였고, detector를 exit 근처 `z=5.2 um`에 놓으면 16.98%, `z=14 um`에서는 47.34%였다. 크기만 줄이고 detector geometry를 그대로 두는 것은 오히려 error를 키웠으므로 size 하나만으로 유효성을 말할 수 없다.

first Rytov는 `U=U0 exp(psi)`를 scalar Helmholtz equation에 대입했을 때 생기는 `grad(psi) dot grad(psi)` 항을 버린다. 큰 누적 phase를 exponential로 보존할 수는 있지만, spheroid 두께 변화와 경계 refraction·diffraction이 만드는 큰 phase gradient까지 정확히 보존하는 것은 아니다. exit 부근에서 만든 Rytov field를 homogeneous water에서 선형 전파한 결과와 각 detector plane에서 Rytov 식을 새로 평가한 결과가 달라졌고, 그 차이는 `z=7 um`에서 6.53%, `z=14 um`에서 26.93%였다.

따라서 이 biological baseline에서는 Born이 크게 실패했고, corrected first Rytov도 vector `Ex`에 대해 약 25%의 complex scattered-field error를 남겼다. 이 시점의 분석만으로는 그 error를 scalar-versus-vector Maxwell reduction과 first-Rytov truncation으로 분리할 수 없었지만, 뒤의 experimental scalar boundary-value reference가 이를 분리했다.

이 분석은 최종 CSV 기준 15개 row를 완료했다. 로그에는 동일한 15개 set을 두 번 실행한 기록이 있으며 서로 다른 30개 case가 아니다. 일곱 row는 voxel과 continuous control을 함께 계산했고, 여덟 row는 larger-FOV·distance·contrast 등의 continuous control만 계산했다. MATLAB Code Analyzer는 변경·신규 MATLAB file 세 개에서 issue 0을 보고했지만, 일반 `git diff --check`는 untracked `forward/`를 검사하지 않으므로 이 tree의 whitespace 검증 증거로 세지 않는다.

이전의 [compare_with_spheroid.m](../forward/compare_with_spheroid.m)은 `x<=20`의 72-case sweep를 정의하지만 실행 결과가 없고, 새 분석은 그 legacy [born_rytov_spheroid.m](../forward/born_rytov_spheroid.m)을 호출하지 않는다. 두 작업을 혼동하면 안 된다.

### 8.5 experimental scalar boundary-value reference와 오차 분리

18:56의 report 뒤 [forward/exp/](../forward/exp/)에서 real prolate, axial incidence, scalar `m=0`에 한정한 scalar Helmholtz boundary-value reference를 별도로 조립했다. 이 reference는 scalar boundary condition과 coefficient solve는 vector solver와 독립적이지만, [sph_eigen.m](sph_eigen.m), [sph_angular.m](sph_angular.m), radial function과 outgoing ODE는 공유하므로 완전히 독립된 numerical implementation은 아니다. [scalar_prolate.m](../forward/exp/scalar_prolate.m)은 surface에서 scalar field `U`와 normal derivative를 연속시켜 외부 outgoing coefficient와 내부 regular coefficient를 푼다. [scalar_prolate_eval.m](../forward/exp/scalar_prolate_eval.m)은 total·scattered field와 analytic Cartesian gradient를 평가하고, [scalar_far_amplitude.m](../forward/exp/scalar_far_amplitude.m)은 큰 공통 phase를 직접 계산하지 않고 outgoing far-amplitude limit를 구한다.

이 solver는 public general-purpose API가 아니라 원인 분리를 위한 experimental reference다. 지원 범위는 real refractive index, prolate, `+z` axial incidence, axisymmetric scalar field뿐이며 oblate, side incidence, arbitrary polarization, absorption은 구현하지 않았다.

첫 near-sphere test는 MATLAB `besselh`의 order와 argument array가 자동 outer expansion되지 않아 실패했다. explicit size expansion으로 고친 뒤 다음 검사를 통과했다.

- zero contrast에서 incident plane wave와 analytic gradient 복원
- aspect ratio가 1에서 `1e-6` 벗어난 prolate와 independent scalar spherical partial-wave solution의 상대차 `1.12e-6`
- independent quadrature grid에서 구한 discrete relative 2-norm boundary residual: near-sphere `1.54e-15`, biological baseline `5.47e-15`
- analytic field gradient와 Cartesian finite difference의 일치
- baseline `L=112`와 `L=132` far amplitude 차이 `1.21e-12`
- analytic far limit와 반경 `1e6 a` evaluation의 상대차 `7.51e-6`
- weak-potential exact derivative와 analytic Born coefficient의 상대차 `5.10e-5`
- exact second coefficient를 추출하는 contrast step refinement 차이 `8.36e-5`

동일한 biological baseline과 ideal circular pupil `NA=0.10`에서 converged scalar-reference scattered field와 full-vector `Ex`의 relative L2 gap은 `0.252653%`였다. first Rytov와 scalar reference의 error는 `z=7 um`에서 `24.174%`이고, 같은 Rytov와 vector `Ex`의 error는 `24.302%`였다. 따라서 이 특정 조건에서는 약 24% residual의 대부분이 vector-to-scalar reduction이 아니라 first-Rytov truncation에서 나온다고 분리할 수 있다.

scattering potential `n_p^2-n_m^2`를 baseline의 비율 `t`로 scale하고 scalar boundary-value reference와 비교한 결과는 다음과 같다. `R2`는 symmetric weak-potential reference solve에서 second Rytov coefficient를 추출한 diagnostic이며, 현재의 일반 forward model 출력이 아니다.

| `t` | `z` (um) | Born error | First Rytov error | Second-order diagnostic error | `norm(t^2 psi2)/norm(t psi1)` |
|---:|---:|---:|---:|---:|---:|
| 0.10 | 5.2 | 16.489% | 1.856% | 0.028% | 0.0317 |
| 0.10 | 7.0 | 16.489% | 2.516% | 0.046% | 0.0421 |
| 0.10 | 14.0 | 16.489% | 4.700% | 0.275% | 0.0731 |
| 0.25 | 5.2 | 41.372% | 4.505% | 0.162% | 0.0792 |
| 0.25 | 7.0 | 41.372% | 6.153% | 0.265% | 0.105 |
| 0.25 | 14.0 | 41.372% | 11.396% | 2.259% | 0.183 |
| 0.50 | 5.2 | 84.860% | 8.663% | 0.707% | 0.158 |
| 0.50 | 7.0 | 84.860% | 12.005% | 1.233% | 0.210 |
| 0.50 | 14.0 | 84.860% | 21.917% | 19.658% | 0.366 |
| 1.00 | 5.2 | 197.233% | 16.873% | 4.393% | 0.317 |
| 1.00 | 7.0 | 197.233% | 24.174% | 13.381% | 0.421 |
| 1.00 | 14.0 | 197.233% | 47.214% | 2,704.996% | 0.731 |

`N=512`와 `N=1024`에서 first-와 second-order error metric의 최대 차이는 각각 `1.28e-9`, `8.86e-8`이었다. second-order term은 약하거나 가까운 조건에서 크게 개선되지만, baseline을 멀리 전파한 `z=14 um`에서는 asymptotic expansion이 무너져 오히려 발산한다. 고차항 하나를 더하는 것이 전 범위의 안정적인 해법이라는 뜻은 아니다.

[diagnose_rytov_term.m](../forward/exp/diagnose_rytov_term.m)은 scalar-reference total field에서 `psi=log(U/U0)`의 생략항 `Q=grad(psi) dot grad(psi)`를 직접 평가했다. baseline `t=1`에서 내부 core의 weighted RMS `|Q|/|f|`는 grid spacing 0.05와 0.025 um에서 `0.2410`, `0.2394`, rim에서는 `0.2863`, `0.2858`이었다. transverse complex-term norm ratio는 core `0.9063–0.9146`, rim `0.9902–0.9903`이었다. `t=0.1`에서는 fine-grid `|Q|/|f|`가 core `0.0149`, rim `0.0330`으로 줄었다. 통계는 log derivative가 불안정한 `|U|<=0.1` 표본을 제외하므로 field zero 주변의 global bound가 아니다. 이 결과는 baseline first Rytov가 버린 quadratic gradient term이 작지 않고, 특히 spheroid의 transverse geometry가 중요하다는 직접적인 진단이다. transverse norm ratio는 complex term의 비율이지 energy fraction이 아니다.

같은 10 um axial thickness, index, wavelength, `z=7 um`, `NA=0.10`, phase delay `4.0886 rad`를 유지한 geometry control에서 prolate first-Rytov scattered-field error는 24.1742%, independent scalar spherical partial-wave reference를 쓴 sphere는 12.7176%였다. uniform slab control은 scattered-field error 2.9774%, total-field error 5.3004%였다. slab 값은 단일 전파 mode의 pointwise ratio라 pupil L2인 prolate·sphere 값과 정량적으로 같은 metric은 아니지만, centerline phase 하나가 Rytov accuracy를 결정하지 않는다는 정성적 control이다.

scalar reference의 pre-objective total field에는 `rho=0.1598748398 um`, `z=6.996322478 um`에서 zero가 발견되었다. `L=112`에서 `|U|=2.41e-14`, `L=132` 재평가에서 `8.64e-12`였고 axisymmetry 때문에 ring zero에 해당한다. 이 점에서는 `psi=log(U/U0)`가 singular하므로 detector plane 주변에서 uniformly small log-gradient를 가정할 수 없다. 이는 pupil-filtered detector field에도 같은 zero가 있다는 주장은 아니며, 이 한 ring이 24.174% error 전부를 설명한다는 뜻도 아니다.

보존된 [diagnosis.log](../forward/exp/diagnosis.log)에는 수치 진단과 함께 실패 이력도 남아 있다. 초기 `besselh` broadcasting failure는 수정되어 test가 통과했고, 처음의 near-sphere focal-coordinate control은 outgoing radial start-refinement 기준을 통과하지 못했다. 최종 geometry control은 그 불안정한 limit 대신 independent scalar spherical partial-wave series를 사용해 완료했다. 그림은 생성하지 않았고, 최종 numerical rerun은 [verification.log](../forward/exp/verification.log)의 `ALL EXPERIMENT CHECKS COMPLETED`로 끝났다.

## 9. 무엇이 성립했고 무엇이 아직 성립하지 않았는가

### 9.1 현재 증거가 지지하는 결론

1. 균질·등방성 prolate와 oblate spheroid의 plane-wave scattering을 Maxwell full-vector separation-of-variables로 계산하는 MATLAB reference solver가 구현되었다.
2. sphere limit, 좌표·특수함수 identity, vector curl/divergence, boundary continuity, Rayleigh limit, optical theorem, lossy energy balance, Barton benchmark를 포함한 검사들을 통과했다.
3. 지정한 biological baseline의 좁은 forward gnomonic patch에서는 산란장의 cross-polarized aggregate norm이 1.995%–3.471%이고, 넓은 patch에서는 12.633%–20.325%다.
4. 따라서 이 biological baseline의 sampled scattered field는 좁은 patch에서 fixed-polarization scalar reduction을 검토할 근거를 주지만, 검사한 wide forward patch에서는 그 운영 기준을 통과하지 못한다.
5. 같은 baseline의 `NA=0.10` detector pupil 안에서 Born의 complex scattered-field error는 197.56%, corrected first Rytov는 25.74%였다. Born은 이 조건과 metric에서 부정확했고, first Rytov도 상당한 residual을 남긴다.
6. relative index contrast `0.001` control에서는 corrected voxel Rytov 1.21%, continuous-shape Rytov 1.04%였다. 근사가 contrast에 민감하다는 것은 확인했지만, 이 두 점만으로 허용 범위의 경계를 정할 수는 없다.
7. experimental scalar prolate boundary-value reference는 baseline에서 vector `Ex`와 0.252653% 차이였고, first Rytov는 이 scalar reference와 24.174% 차이였다. 이 제한된 case에서는 scalar polarization reduction보다 first-Rytov truncation이 지배적인 오차원이다.

### 9.2 아직 성립하지 않은 주장

1. Born과 first Rytov의 복소 field error는 일부 조건에서 측정했지만, wavelength·contrast·size·aspect ratio·orientation·NA 전반의 유효 범위는 정해지지 않았다. low polarization leakage는 scalar model accuracy의 충분조건이 아니다.
2. 현재 biological geometry의 phase-delay heuristic은 baseline에서 `4.089 rad`, half-size control에서 `2.044 rad`였다. 둘 다 plan에서 고려한 `pi/2` Born heuristic보다 크며, 실제 baseline Born error 197.56%와도 일관된다.
3. [tests/test_scalar_approximation.m](tests/test_scalar_approximation.m)은 aggregate Frobenius norm을 사용하므로 모든 sampled direction에서 co-polarized component가 지배적이라는 pointwise 보장을 주지 않는다.
4. polarization test와 Born/Rytov comparison은 scattered field를 평가한다. total 또는 transmitted detector field의 scalar reduction은 검증하지 않았으며, total-field 편광은 강한 incident field 때문에 지배적으로 보일 수 있어 다른 질문이 된다.
5. wide-patch leakage에는 `E_sca dot rhat=0`인 far-field transversality가 강제하는 geometric projection이 포함된다. 따라서 이 수치만으로 scalar Helmholtz·Born·Rytov model의 amplitude 실패를 판정할 수 없다.
6. polarization test의 gnomonic square sampling은 solid-angle weighted NA integral이 아니다. Born/Rytov comparison은 이상적인 circular pupil `NA=0.10`을 사용하지만 실제 objective transfer function, polarization analyzer, detector response는 포함하지 않는다.
7. Born/Rytov 분석에는 weak contrast, half size, 633 nm, detector distance와 numerical-resolution control이 있지만 systematic parameter sweep는 아니다. oblate, side incidence, y/z polarization도 complex-field comparison에서 다루지 않았다.
8. 세포는 실제 세포가 아니라 homogeneous, isotropic, nonabsorbing single-index spheroid다. membrane, nucleus, organelles, spatial heterogeneity, dispersion, absorption, culture-medium composition을 포함하지 않는다.
9. plane-wave illumination만 지원한다. focused beam과 illumination NA는 구현되지 않았다.
10. spheroid에 대한 독립 문헌 비교는 한 aspect ratio와 incidence에서 얻은 Barton의 네 far-field intensity maximum뿐이다. complex phase, 전체 angular pattern, near field, internal field에 대한 독립 spheroid benchmark는 실행하지 않았으며, runtime residual·Maxwell identity·energy·Rayleigh·cutoff 검사는 강한 consistency evidence이지만 이를 대체하지 않는다.
11. experimental scalar reference로 baseline의 오차 분리는 수행했지만, 그 solver는 real prolate axial case에 한정되고 fixed `L`을 사용한다. scalar boundary assembly와 sphere check는 별도지만 vector solver의 eigen·angular·radial·outgoing-ODE 구현을 공유하므로 correlated special-function error 가능성이 남고, production solver의 adaptive convergence contract를 갖추지 않았다.
12. corrected voxel과 continuous-shape Rytov 계산은 forward-propagating mode만 포함하며 evanescent와 exactly grazing mode를 제외한다. 따라서 `z=5.2 um` 결과는 surface에서 0.2 um 떨어진 propagating-field control이지 evanescent near-field validation이 아니다.
13. oblate, side incidence, 다른 polarization, absorbing particle에서는 scalar-versus-vector와 approximation error가 아직 분리되지 않았다.
14. 이전 72-case `x<=20` comparison sweep와 전체 537-case vector-solver stress campaign은 실행되지 않았다. 완전한 11-suite one-tree run도 아직 없다.
15. corrected Rytov operator 뒤 기존 inverse test와 default demo는 갱신되지 않았다. 기존 `dx=0.2 um` 설정은 현재 `born_rytov_voxel:RytovNyquistViolation`에서 멈추며, pupil과 logarithm이 commute한다고 가정한 과거 machine-precision equality도 더 이상 유효하지 않다.
16. high-order complex oblate eigenproblem 일부는 vector-solver parameter envelope 안에서도 representation-sensitive하여 명시적으로 error를 낸다.
17. `xi0=0`인 limiting oblate geometry는 outgoing log-ODE normalization 때문에 지원하지 않는다.

## 10. 다음 단계의 우선순위

full-vector co-polarized field와 Born/Rytov complex scattered field의 직접 비교, 그리고 한 prolate baseline에서 experimental scalar boundary-value reference를 이용한 오차 분리는 완료되었다. 다음 핵심은 이 reference를 독립 test가 있는 재현 가능한 solver로 정리하고 oblate·orientation·polarization·absorption까지 일반화한 뒤, scalar error의 경계를 systematic하게 찾는 것이다. 최종 목적이 검출된 `Ex` 자체의 예측이라면 vector-reference error가 end-to-end metric이고, 근사 실패의 원인을 분석할 때 converged scalar boundary-value comparison이 필요하다.

그 다음에는 다음 순서가 타당하다.

1. experimental scalar prolate code를 정식 test 대상에 포함하고, oblate와 non-axial incidence로 확장하기 전에 independent benchmark와 convergence contract를 보강한다.
2. propagation direction을 고정한 채 두 transverse polarization을 비교하고, orientation sweep는 별도 변수로 분리한다.
3. wavelength, contrast, absorption, size, aspect ratio, detector distance, NA를 sweep하여 scalar error map과 목표 error threshold를 만든다.
4. corrected Rytov convention과 sampling limit에 맞게 inverse test와 demo를 갱신한다.
5. 전체 11-suite를 현재 tree에서 한 번 실행하고, 이후 537-case stress campaign을 완료한다.
6. 실험을 목표로 한다면 homogeneous-cell 모델 다음에 layered 또는 heterogeneous numerical reference를 별도 solver로 검증한다.

## 11. 재현 가능한 증거와 실행 주의사항

검증 증거는 저장소 바깥의 임시 디렉터리와 현재 untracked `forward/` tree에 나뉘어 있다.

- `/Users/jhlee/Documents/jhlee_forge/projects/tmp/spheroid-analytic-forward/sdd/`: 단계별 구현 report
- `/Users/jhlee/Documents/jhlee_forge/projects/tmp/spheroid-analytic-forward/runtime/`: red/green 및 focused runtime log
- `/Users/jhlee/Documents/jhlee_forge/projects/tmp/spheroid-analytic-forward/green/run_all_final10.log`: commit `569307b`의 10-suite 기준 로그
- `/Users/jhlee/Documents/jhlee_forge/projects/tmp/spheroid-analytic-forward/green/checkcode_final5.log`: 30개 MATLAB file, Code Analyzer issue 0
- `/Users/jhlee/Documents/jhlee_forge/projects/tmp/spheroid-analytic-forward/green/example_xz_slice_final5.log`: example field 실행 결과
- `/Users/jhlee/Documents/jhlee_forge/projects/tmp/spheroid-analytic-forward/green/stress_envelope_smoke_final.log`: stress report/resume smoke test
- `/Users/jhlee/Documents/jhlee_forge/projects/tmp/spheroid-analytic-forward/biological_scalar.log`: biological scalar-polarization 6-case 결과
- `forward/tests/rytov_analysis.md`, `.csv`, `.log`: corrected detector-NA Rytov와 15-row vector-reference 비교
- `forward/exp/report.md`: experimental scalar-reference 분석의 최종 해석과 재현 절차
- `forward/exp/verification.log`: 최종 실험 rerun; `ALL EXPERIMENT CHECKS COMPLETED`
- `forward/exp/diagnosis.log`: preliminary probe와 수치 진단을 포함한 개발 이력
- `forward/exp/contrast_distance_N512.csv`, `contrast_distance_N1024.csv`: scalar reference에 대한 Born, first-와 second-order Rytov error
- `forward/exp/omitted_term_dx0.1.csv`, `omitted_term_dx0.05.csv`, `omitted_term_dx0.025.csv`: 생략된 `grad(psi) dot grad(psi)`의 grid-refinement 통계
- `forward/exp/geometry_control.csv`: 같은 phase-delay의 prolate와 sphere 비교
- `forward/exp/field_zero.csv`: pre-objective scalar total-field zero의 위치와 refinement 값

최종 command들은 `-batch` 대신 직렬화된 `-nojvm -nodesktop -nosplash -r` 방식으로 실행했고, session-local `restoredefaultpath` 뒤 test는 exit status 0으로 완료되었다.

현재 [tests/run_all.m](tests/run_all.m)은 scalar-polarization test를 포함한 11개 suite를 나열한다. 그러나 보존된 10-suite 기준 로그와 약 2시간 11분의 scalar 로그는 서로 다른 실행이므로, 새 11-suite 전체 로그가 생기기 전에는 "11개 모두 한 번에 통과"라고 기록하지 않는다. `forward/tests/`와 `forward/exp/`의 실험도 이 runner와 별도이며 모두 untracked 상태다.

## 12. 최종 결론

이 작업에서 완성된 것은 특정 균질 spheroid의 full-vector scattering field를 계산하고 실행별 수치 검증을 수행하는 reference solver다. 가장 어려웠던 부분은 공식 자체보다 coordinate singularity, high-order scaled special functions, outgoing radial basis normalization, ill-conditioned modal systems을 floating-point arithmetic에서 fail-closed 방식으로 다루는 일이었다.

현재 scalar 결과는 한 biological baseline에서 원인을 수치적으로 분리했다. real prolate, axial x-polarized plane wave, ideal `NA=0.10`의 complex scattered-field L2 metric에서 scalar boundary-value reference와 vector `Ex`의 차이는 0.252653%였고, Born과 first Rytov의 scalar-reference error는 각각 197.233%, 24.174%였다. 이는 이 한 조건에서는 first-Rytov truncation이 주된 오차원이라는 뜻이며, 일반 유효 범위와 실제 heterogeneous cell·microscope model의 정확도는 아직 검증되지 않았다.

## 13. 2026-09-12 추가 기록: 제안된 oblate 구조의 scalar approximation과 NA

이 절은 위 기록 이후 같은 대화에서 수행한 실험과 설명을 추가한다. 앞 절의 결론과 미완료 항목은 당시의 기록으로 보존하며, 아래 결과가 새로운 oblate 조건의 scalar-versus-vector 비교 범위를 보완한다. 이번 계산은 MATLAB session 안에서만 수행했고, 기존 solver나 `forward/exp/`의 scalar prolate API를 수정하지 않았다. 아래 수치는 대화에 남은 완료된 console 출력을 정리한 것이며, 이 문서 추가 시점에 재실행한 결과는 아니다.

### 13.1 변경한 물리 조건과 조사 범위

| 항목 | 기존 prolate baseline | 이번 oblate 실험 |
|---|---:|---:|
| x, y 방향 전체 길이 | 각각 5 um | 각각 10 um |
| z 방향 전체 두께 | 10 um | 6 um |
| Semi-axes `(a,b)` | `(5,2.5) um` | `(3,5) um` |
| 대상 굴절률 `n_p` | 1.370 | 1.365 |
| 배경 굴절률 `n_m` | 1.335381534 | 동일 |
| 진공 파장 | 0.532 um | 동일, 532 nm |
| 관측면 z | 7 um | 5 um |
| 물체의 +z 끝에서 관측면까지 | 2 um | 2 um |

배경 매질 내 파장은 약 `398.388 nm`, 새 물체의 부피는 `314.159 um^3`, 중심선을 지나는 ray의 phase-delay estimate는 `2.09885 rad`다. 평균 핵 부피 약 `300 um^3`는 사용자가 제안 근거로 전달한 값으로, 이 작업에서 원논문을 특정하거나 검증하지 않았다. 축 길이는 논문에서 추출한 값이 아니라 실험용 설정이며, 모델은 기존 배경 안의 균질·등방성·비흡수성 핵 모양 물체로서 세포질이나 층상 경계를 포함하지 않는다.

처음에는 정상입사에서 detector NA만 바꿨다. 이후 사용자가 비정상입사도 확인하도록 요청하여 검출기 광축은 `+z`에 고정하고 입사 방향을 바꿨다. illumination NA는 처음 전달된 `0.9`를 사용자의 정정에 따라 `0.5`로 대체했으며, 최종 범위는 `NA_illum=n_m sin(theta_inc)<=0.5`, 즉 최대 입사각 `21.98877 deg`다. 여기서 illumination NA는 개별 tilted plane wave의 입사각 범위를 뜻하며, focused beam이나 여러 입사각의 동시·비간섭 조명을 구현한 것은 아니다.

### 13.2 x polarization의 정의

[incident_plane_wave.m](incident_plane_wave.m)의 `pol=[1,0]`을 그대로 사용했다. `theta=theta_inc`, `phi=phi_inc`일 때 convention은 다음과 같다.

```text
khat = (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta))
e_TE = (-sin(phi), cos(phi), 0)
e_TM = (cos(theta) cos(phi), cos(theta) sin(phi), -sin(theta))
e0   = cos(phi) e_TM - sin(phi) e_TE
     = Rz(phi) Ry(theta) Rz(-phi) xhat
```

이는 정상입사의 x 방향 전기장을 입사 방향과 함께 추가 roll 없이 회전시키는 정의다. 정상입사에서는 `e0=xhat`, `phi=0`에서는 `e0=cos(theta) xhat-sin(theta) zhat`이며 항상 `e0 dot khat=0`이다. 입사 방향만으로 편광이 유일하게 정해지는 것은 아니므로 이 convention이 필요하다. 정상입사를 고정하고 detector NA를 늘리는 경우에는 입사 전기장을 회전시키는 것이 아니라 더 다양한 산란 방향을 수집하며, 각 방향의 산란 편광은 Maxwell 해가 결정한다.

### 13.3 비교 대상과 오차 정의

새 실험에서는 Born이나 first Rytov가 아니라, 표면에서 `U`와 normal derivative를 연속시킨 scalar Helmholtz boundary-value 해를 reference로 사용했다. oblate와 tilted incidence에 필요한 scalar system을 session 안에서 조립했으며, vector 해와 eigen·angular·radial 특수함수 및 outgoing ODE를 공유하므로 완전히 독립된 구현은 아니다.

같은 ideal circular pupil을 통과한 complex scattered field를 각각 vector `E_s`, scalar `U_s`라 두고 다음 relative L2 metric을 계산했다.

```text
u_co         = E_s dot conj(e0)
epsilon_perp = ||E_s-u_co e0||_2 / ||E_s||_2
epsilon_co   = ||U_s-u_co||_2 / ||u_co||_2
epsilon_full = ||U_s e0-E_s||_2 / ||E_s||_2
epsilon_full^2 = epsilon_perp^2 + (1-epsilon_perp^2) epsilon_co^2
```

`epsilon_perp`는 고정 입사 편광 밖의 성분, `epsilon_co`는 co-polarized complex amplitude의 scalar-model 오차, `epsilon_full`은 scalar field를 `e0` 방향 vector로 복원했을 때의 전체 오차다. intensity, total field, 영상 복원 오차가 아니다.

적분은 propagating angular spectrum과 Parseval 관계로 전체 검출기 평면의 L2 norm을 구했다. `E_s(r) ~ A(khat) exp(i k_m r)/r`인 far amplitude로부터 `S(kx,ky;z)=2 pi i A(khat) exp(i kz z)/kz`를 사용했으며, radial NA 적분의 공통 상수를 제외한 weight는 `NA/(n_m^2-NA^2)`다. 따라서 이전 gnomonic square의 unweighted norm이나 solid-angle power 적분과 다른 metric이다. 공통 propagation phase는 이 norm에서 상쇄되므로, `z=5 um` 설정은 유지했지만 이 결과를 유한 camera window나 evanescent near field의 거리 의존성을 검증한 것으로 해석하면 안 된다.

### 13.4 정상입사 pilot

정상입사 vector 해는 public `spheroid_solve`의 `tol=1e-6` runtime validation을 통과했고 `info.validated=true`, `L=109`, `M=11`이었다. scalar는 `L=112`와 `132`를 비교했으며 detector NA 간격은 `0.0005`, 최대값은 `1.3`이었다. 아래 값의 단위는 모두 %다.

| Detector NA | `epsilon_perp` | `epsilon_co` | `epsilon_full` |
|---:|---:|---:|---:|
| 0.028 | 0.996217 | 0.060604 | 0.998059 |
| 0.05 | 1.544719 | 0.058108 | 1.545811 |
| 0.10 | 1.902210 | 0.075638 | 1.903713 |
| 0.30 | 2.478660 | 0.174686 | 2.484804 |
| 0.50 | 2.715368 | 0.303327 | 2.732245 |
| 1.00 | 3.023760 | 0.726344 | 3.109697 |
| 1.20 | 3.154261 | 1.079322 | 3.333637 |
| 1.30 | 3.263472 | 1.468017 | 3.578132 |

정상입사에서 전체 vector 오차의 첫 1% crossing은 detector NA 약 `0.0281`로, 실제 grid bracket은 `[0.0280,0.0285]`였다. co-polarized scalar-model 오차만의 첫 1% crossing은 약 `1.1755`로 `[1.1755,1.1760]` 안에 있었다. 따라서 co-polarized amplitude가 정확하다는 주장과 전체 vector field를 1% 이내로 축약할 수 있다는 주장은 구분해야 한다.

### 13.5 illumination NA 0.5까지 확장한 결과

입사 polar angle은 `NA_illum=0:0.05:0.5`의 11개와 정확히 `15 deg`인 추가점, 총 12개를 사용했다. 각 polar angle에서 incidence plane의 TE/TM 해를 구한 뒤 회전대칭과 선형성을 이용해 위 x-polarization을 구성했다. `phi_inc=0,45,90 deg` 혼합을 직접 평가했고, 원형 pupil 적분에서 TE/TM cross term이 사라지는 reflection symmetry를 수치 확인하여 모든 입사 azimuth의 `epsilon_full`이 두 endpoint 사이에 놓임을 확인했다. polar angle은 이산 표본이므로 연속된 모든 입사각에 대한 증명은 아니다.

계산은 signed m마다 동일한 Galerkin boundary matrix를 한 번 조립하고 여러 입사각·편광의 right-hand side를 함께 푸는 batch 방식으로 수행했다. coarse cutoff는 `(L,M)=(109,55)`, fine은 `(114,59)`였고, 입사장 azimuth quadrature는 256점, 관측 azimuth는 128점, detector NA는 `0:0.001:1.3`이었다. 아래는 조사한 입사 조건에 대한 각 metric의 최대값이며 모두 %다. 열별 최대값은 서로 다른 입사 조건에서 나올 수 있으므로 최대값끼리 위 오차 분해식을 적용하면 안 된다.

| Detector NA | 최대 `epsilon_perp` | 최대 `epsilon_co` | 최대 `epsilon_full` |
|---:|---:|---:|---:|
| 0.10 | 35.979103 | 13.078857 | 37.992225 |
| 0.30 | 25.052724 | 7.395468 | 26.055694 |
| 0.50 | 2.933348 | 0.527766 | 2.980407 |
| 0.60 | 2.795397 | 0.417567 | 2.820432 |
| 0.70 | 2.862000 | 0.469306 | 2.897349 |
| 0.90 | 3.004018 | 0.623498 | 3.064197 |
| 1.10 | 3.171961 | 0.871303 | 3.273002 |
| 1.30 | 3.448761 | 1.478862 | 3.696564 |

detector NA `0.5–1.3`의 801개 표본 모두에서 최악의 `epsilon_full`은 5% 이내였고, 구간 전체 최대는 약 `3.697%`였다. 따라서 이 모델·파장·metric과 조사한 illumination 범위에서는 5% 및 기존 10% 허용 기준을 통과하지만, 1% 기준은 만족하지 않는다. 이는 NA가 큰 pupil에 포함된 모든 개별 산란 방향이 정확하다는 뜻이 아니라, 전체 pupil에 걸친 field norm의 결과다.

작은 detector NA가 비정상입사에서도 항상 유리한 것은 아니었다. 입사각이 커지면 강한 forward-scattering lobe가 고정된 검출기 광축에서 벗어나므로, 좁은 pupil은 상대적으로 큰 산란각의 약한 field를 수집하여 relative error가 커질 수 있다. 특히 `NA_det<NA_illum`이면 직접 입사광 방향도 pupil 밖에 있으므로 일반적인 transmitted-reference 측정과 구분해야 한다. 정상입사의 작은-NA 결론을 illumination NA 0.5 전체에 그대로 적용할 수 없다.

### 13.6 수치 검증과 실행 중 발생한 문제

정상입사 scalar pilot은 다음 검사를 통과했다.

- scalar `L=112`와 `132`의 far-amplitude 상대차: `3.23e-13`
- scalar boundary residual: 각각 `3.35e-15`, `4.40e-15`
- zero-contrast scattered far-amplitude norm: 약 `1.87e-15`
- weak-potential symmetric derivative와 analytic Born coefficient의 상대차: `4.92e-7`
- 거의 구인 oblate와 independent scalar spherical partial-wave 해의 상대차: `1.54e-6`
- vector far amplitude의 두 Richardson extrapolation 간 상대차: `6.64e-8`

비정상입사 batch의 완료된 검증은 다음과 같다.

| 검증 | 결과 |
|---|---:|
| coarse/fine cutoff vector field 상대차 | `1.82e-12` |
| coarse/fine cutoff scalar field 상대차 | `2.08e-12` |
| coarse/fine cutoff 최대 metric 차이 | `8.29e-13` |
| 별도 boundary grid의 최대 vector residual | `7.02e-11` |
| 별도 boundary grid의 최대 scalar residual | `1.68e-11` |
| analytic vector far limit와 public evaluator의 TM 상대차 | `2.95e-8` |
| 같은 비교의 TE 상대차 | `3.02e-8` |
| NA 간격 0.001/0.002의 최대 metric 차이, NA >= 0.5 | `8.34e-6` |
| 원형 pupil 적분의 normalized TE/TM cross term | `9.79e-15` |

별도 boundary grid는 `248 theta x 256 phi`였고, 모든 12개 입사각과 두 편광에서 vector의 tangential E/H·normal `n^2 E`, scalar의 U·normal derivative 연속성을 확인했다. far-limit 비교는 최대 입사각에서 다섯 관측 방향의 public `spheroid_eval` 결과를 큰 반경과 두 배 반경으로 extrapolate하여 수행했다. NA quadrature의 `8.34e-6`은 무차원 error ratio 차이로 약 `0.000834 percentage points`이며, cutoff 간 작은 차이 자체가 절대 정확도 보증은 아니다. 거의 0인 high-m right-hand side로 나눈 sector별 상대 residual은 커질 수 있어, 위 표에는 재구성된 전체장의 독립 grid residual을 기록했다.

public `spheroid_solve`로 별도 시작한 `15 deg`와 최대 입사각의 TM 계산은 validation 단계가 길어져 batch 검증 완료 후 수동 중단했다. 따라서 두 public solve의 `info.validated=true`를 얻었다고 기록하지 않는다. batch 검증은 위 별도 검사들로 완료했으며, 최종 console marker는 `ALL_BATCH_EXPERIMENT_CHECKS_PASS`였다.

MATLAB R2024a Update 9의 no-JVM session에서 실행했다. session 내부 default path 복구로 진행했으며 사용자 전역 설정은 수정하지 않았다. 성공 증거로는 정상 실행과 후속 assertion 결과만 사용했다. 소유한 세 MATLAB session은 작업 후 모두 종료되었다.

### 13.7 결론과 보존 상태

이번 결과는 제안된 oblate 구조에서 scalar Helmholtz 해와 Maxwell vector 해의 차이를 정상입사 및 illumination NA 0.5까지의 이산 polar-angle sweep으로 평가한 것이다. scalar approximation의 성립 여부는 허용 오차와 함께 명시해야 하며, ideal detector NA `0.5–1.3`에서는 조사 범위의 전체 complex scattered-field 오차가 5% 이내였으나 1% 이내는 아니었다. 새로운 oblate 조건의 Born·Rytov truncation error, 실제 objective의 polarization transport·aberration·analyzer, finite camera window, evanescent field, 연속 polar-angle bound, 다른 파장·굴절률·형상은 검증하지 않았다.

이 추가 기록 이전에는 이번 oblate 실험의 `.m`, CSV, MAT, 별도 log 파일을 저장하지 않았고, session 종료로 계산 배열도 남아 있지 않다. 따라서 이 절은 대화의 console 증거를 보존하는 문서이지 실행 가능한 재현 artifact가 아니며, 기존 `forward/exp/scalar_prolate.m`이 oblate·tilted incidence를 지원하게 되었다는 뜻도 아니다. 후속 재현이 필요하면 session에서 사용한 batch 조립·검증 절차를 최소 script로 정리하고 다시 실행하여 입력·결과·로그를 함께 저장해야 한다.
