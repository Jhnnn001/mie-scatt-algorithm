# Δn × 0.01 고 NA inverse reconstruction 보고서

상태: 2026-09-13, 예정한 11개 run·17개 복원과 sampling·민감도 검사를
완료했다. 미수렴 run도 결과로 보존했으며, 완료를 수치 수렴과 혼동하지 않는다.

## 요약 및 결론

선택한 Δn × 0.01 시스템에서 낮은 scalar forward 오차는 확인했지만,
현재 세 방법으로 정확한 3D 굴절률 복원이 확보되지는 않았다.
아래는 같은 padding 8 연산자를 사용한 비교다. GP는 기본 100회 cap,
TV는 물리적으로 penalty를 맞추고 tolerance 10⁻⁵로 검사한 결과다.

| 방법 | RI 오차 | 물리적 산란장 오차 | 수치 상태 |
|---|---:|---:|---|
| Direct Fourier | 46.9056% | 17.6344% | 직접 변환 |
| 비음수 GP | 36.3767% | 0.2263% | 100회 cap, 미수렴 |
| Lim TV, strict | 28.9450% | 0.9293% | 856 inner / 1 outer, 수렴 |

TV가 이 비교에서 가장 낮은 RI 오차를 보였지만, 굴절률 과소평가와
축방향 경계 artifact가 남았다. Padding 4의 TV 오차 19.5501%는 더
좋은 단일 수치일 뿐, refinement에 안정적인 검증 성능으로 선택할 수 없다.
또한 padding 4에서 GP를 300회·더 엄격한 CGLS로 실행한 결과는
85.0047% RI 오차와 미수렴으로 악화되었다.
따라서 단순한 반복 증가나 Δn 감소만으로 현재 inverse 문제가 해결된다는
결론은 틀리며, missing cone, 데이터·모델 불일치, 수치 projection과
정규화·경계 조건을 구분해서 다뤄야 한다.

## 1. 대상 시스템

| 항목 | 값 |
|---|---|
| 진공 파장 λ₀ | 0.532 μm |
| 물 굴절률 nₘ | 1.335381534 |
| 입자 굴절률 nₚ | 1.33567771866 + 0i |
| Δn = nₚ − nₘ | 0.00029618466 |
| 형상 | oblate, a = 3 μm(z), b = 5 μm(x, y) |
| 전체 크기 | 10 × 10 × 6 μm |
| 조명 / 검출 NA | 0.5 / 1.0 |
| 매질 내 최대 각도 | 21.989° / 48.491° |
| 조명 | 81방향, kᵢ,z > 0 |
| 검출면 | z = +5 μm, 고정 평면 |
| 검출 격자 | 512 × 512, 간격 0.1 μm, 폭 51.2 μm |
| 복원 격자 | 128 × 128 × 80, 간격 0.1 μm |
| 잡음 | 추가하지 않음 |

따라서 횡방향 복원은 128²이며, 2차원 inverse가 아니라 128 × 128 × 80의
3차원 굴절률 복원이다. 검출 격자와 복원 격자를 독립적으로 다룬다.
빛은 −z 쪽에서 입사하며, 비스듬한 조명에서도 검출면은 회전하지 않는다.
정답 spheroid의 모양, 균질성, 실제 support 및 굴절률 상한은 inverse에
주지 않았고, 실제 정답은 오차 평가와 별도로 표시한 matched control에만 썼다.
유한 계산 상자와 n ≥ nₘ 조건은 여전히 복원 모델의 가정이다.
선택한 nₚ는 Δn를 인위적으로 0.01배 낮춘 검증용 값이며, 이 낮은 값을
일반적인 실제 생세포 굴절률이라고 해석하지 않는다.

## 2. exact 데이터와 scalar 모델의 범위

완료된 [forward contrast 실험](../../spheroid-analytic-forward/2/readme.md)의
`case_3_reference.mat`를 사용했으며, 별도 voxel 모델로 exact 관측값을
합성하지 않았다. 기존 forward helper의 편광과 Maxwell co-polarized
projection을 그대로 사용하고, angular spectrum에 실제 검출 pupil을
적용한 뒤 total field U를 구성했다.

이 contrast에서 기존 forward의 각도 가중 평균 오차는 scalar Rytov 대
Maxwell co-polarized field가 0.658709%, scalar exact 대 Maxwell co가
0.647658%이며, scalar Rytov 대 scalar exact는 0.090536%이다.
이 세 값은 본 보고서의 역산 잔차와 정규화가 다르므로 서로 같은 수치로
취급하면 안 된다. exact reference의 cutoff 오차는 vector 4.68 × 10⁻¹¹,
scalar 5.04 × 10⁻¹¹이며, public far-field 교차검사 차이는 4.32 × 10⁻⁸이다.

여기서 확인한 것은 기존 편광·co-projection에 대한 scalar 모델의 정확도다.
임의의 편광에 대한 보편적 정확도나 실제 고 NA objective의 vector
polarization transfer까지 검증한 것은 아니다. 물리적 입사장은 항상
E₀·kᵢ = 0을 만족해야 하며, 비스듬한 입사에서 고정 Cartesian x 벡터를
그대로 쓰면 이 조건을 어길 수 있다.

## 3. forward와 inverse의 수학적 정의

χ = n²/nₘ² − 1, f = kₘ²χ, n = nₘ√(1 + χ).

χ̂(q) = ∫ χ(r) exp(−iq·r) d³r, q = kₛ − kᵢ.

Rytov 전처리는 실제 관측 total field에 대해 L = u₀ log(U/u₀)를
계산하고, 2차원 Fourier transform을 통해 다음 샘플을 만든다.

g(q) = L̂(kₛ,⊥)(−2ikₛ,z) exp(−ikₛ,z z_det)/kₘ².

실제 pupil 연산을 P라 하면 forward는 먼저 pupil 이전의 B와 ψ = B/u₀를
계산한 뒤 P[u₀ expm1(ψ)]를 사용한다. 일반적으로 log와 P는 교환되지
않으므로, pupil을 통과한 U에서 얻은 g = Aχ는 근사식이지 항등식이 아니다.
검출 NA가 pupil에서 생긴다는 설명은 맞지만, pupil을 ψ에 미리 적용하는
것과 복소장에 적용하는 것은 서로 다른 forward 모델이다.
현재 코드는 복소장에 pupil을 적용하는 기존 수정 forward와 일치한다.

이산 A는 원래 Ewald 좌표에서 trilinear interpolation하는 padded Fourier
연산자이며, adjoint는 그 정확한 복소 전치다. 관측 interpolation에 필요
없는 주파수만 제외하는 분리 DFT/FFT 계산은 전체 3차원 FFT와 일치한다.
조명 가중치는 기존 실험의 pupil 면적 quadrature를 쓰며 합은 1이다.
전처리의 두 phase-unwrapping 순서 차이는 0, 최대 경계 위상은
1.595 × 10⁻⁵ rad였다.

### 비교한 세 가지 방법

1. **Direct Fourier:** native 복원 격자에 trilinear deposition하고,
   density normalization·Hermitian completion·zero fill 뒤 inverse FFT한다.
   이 별도 native gridding은 A의 정확한 pseudoinverse가 아니다.
2. **비음수 GP:** CGLS로 minimum-norm least-squares data correction을
   근사한 뒤 χ ≥ 0으로 투영한다. Cartesian bin 교체 대신 실제 비균일
   Ewald 연산자를 사용하며, CGLS normal residual은 해 오차의 상한이 아니다.
   데이터가 모델과 불일치할 때 이 두 투영의 반복은 비음수 least-squares
   최소화와 동일하지 않으며, data 잔차가 단조 감소한다는 보장도 없다.
3. **Lim TV:** 아래 inner 최소화와 outer data reinjection을 사용한다.

χₖ₊₁ = argmin(χ ≥ 0) { α TV(χ) + ½‖Aχ − gₖ‖₂² }.

gₖ₊₁ = gₖ + y − Aχₖ₊₁, g₀ = y.

TV(χ) = Σⱼ √[(Dₓχ)ⱼ² + (Dᵧχ)ⱼ² + (D_zχ)ⱼ²].

Lim 논문의 식 (14), (15)와 Appendix A의 gradient/positivity splitting을
사용한다. 현재 A* A는 diagonal Fourier mask가 아니므로 PCG로 풀고,
논문의 고정 iteration 수 대신 residual과 KKT stationarity를 검사한다.
차분 경계는 periodic이다. 첫 outer에서 종료하면 fixed-data penalized
TV 해를 반환하며, 두 번째 outer부터 data reinjection이 반영된다.

Aχ = c√W χ̂_h, y = c√W g로 둘 때 α = λ c² dx²를 사용한다.
λ = alpha_relative × max|χ_direct|이며, true shape나 true Δn로 λ를 맞추지
않았다. dx² TV(χ)가 물리적 TV 적분을 근사하므로, 이 변환은 padding과
voxel 간격을 바꿀 때 data/TV의 상대적 크기를 유지한다.
이 λ의 숫자는 μm 단위와 현재 이산 주파수 합에 종속되며 보편 상수가 아니다.
동일한 splitting penalty라도 유한 iteration에서의 진행률까지 불변인
것은 아니므로, 수치 수렴과 iteration 민감도를 별도로 확인해야 한다.

## 4. 오차의 정의와 해석

ε_n = ‖n_rec − n_true‖₂ / ‖n_true − nₘ‖₂.

분모에 큰 배경 굴절률을 넣으면 아주 작은 Δn의 실패를 가릴 수 있으므로,
전체 계산 부피에 대한 contrast-normalized RI 오차를 사용한다.
물체 내부 오차와 배경 오차는 같은 분모로 나누며,
ε_n² = ε_inside² + ε_background²를 만족한다.

ε_data = ‖Aχ_rec − y‖₂ / ‖y‖₂.

ε_field = ‖√W(U_s,pred − U_s,exact)‖₂ / ‖√W U_s,exact‖₂.

마지막 값은 total field가 아닌 **산란장**에 대한 오차이며, reconstructed
χ를 다시 물리적 Rytov forward에 넣어 Maxwell co 산란장과 비교한다.
matched control의 data 잔차는 합성 선형 데이터에 대한 값이고, 그 control의
physical-field 오차는 여전히 독립 exact 데이터와 비교한 값이다.
축방향 길이는 중심선에서 true Δn의 절반을 넘는 첫·마지막 voxel의 span으로
정의한다. voxel truth는 5.9 μm이고 연속 spheroid의 실제 직경은 6 μm다.
`axial_span_censored`가 참이면 half-contrast 영역이 계산 상자 경계에
닿은 것이므로, 측정된 span을 완전히 분해된 물체 두께로 해석하지 않는다.

작은 ε_data가 작은 ε_n을 보장하지 않는다. missing cone에서는 일부
구조 변화가 관측장을 매우 조금만 바꾸므로 A의 작은 singular value에
대응하는 방향이 불안정하다. 따라서 모델 오차가 작은 경우에도 역산에서
증폭될 수 있으며, 유한 반복과 정규화 오차를 이 효과와 구분해야 한다.
Δn를 낮추면 forward 근사의 정확도는 좋아질 수 있지만 관측각과
missing cone은 변하지 않는다. 선형 χ 모델에서는 A(sχ) = sAχ이므로,
신호 크기만 줄여도 상대적인 복원 오차가 자동으로 줄어드는 것은 아니다.

## 5. 수치 sampling 점검

[operator_convergence.csv](operator_convergence.csv)에서 연속 spheroid의
해석적 χ̂와 voxel χ̂를 비교했고, 같은 χ̂와 실제 Rytov 데이터도 비교했다.
연속 정답의 관측 데이터 mismatch는 0.569585%이다.

| Padding | dx (μm) | Voxel 대 연속 χ̂ | Voxel 대 관측 g |
|---:|---:|---:|---:|
| 2 | 0.10 | 6.330283% | 6.229133% |
| 3 | 0.10 | 3.139894% | 3.063488% |
| 4 | 0.10 | 1.689156% | 1.674612% |
| 6 | 0.10 | 0.952972% | 1.016395% |
| 8 | 0.10 | 0.631425% | 0.792882% |
| 3 | 0.08 | 3.049929% | 2.969663% |
| 6 | 0.08 | 0.818344% | 0.886374% |

0.08 μm에서는 동일한 부피의 160 × 160 × 100 격자를 사용했다.
Padding 2의 오차를 단순히 scalar 근사 실패나 missing cone 때문이라고
해석하면 잘못이다. 이 경우 Fourier interpolation 오차가 지배적이다.
반면 padding을 늘리는 것은 새 관측각을 추가하지 않으므로 missing cone을
없애지 않는다. 위 표는 forward sampling 검사이며 그 자체로 RI 복원의
grid convergence를 증명하지는 않는다.

[detector_convergence.csv](detector_convergence.csv)의 7개 대표 조명에서,
512²·0.1 μm 대비 1024²·0.05 μm의 동일 window/동일 q 데이터 차이는
2.87 × 10⁻¹¹, 1024²·0.1 μm로 window를 두 배 늘린 차이는
1.35 × 10⁻⁶이었다. 따라서 검사한 조건에서는 검출 pitch/window 오차가
Fourier interpolation 오차보다 훨씬 작다.

## 6. 기본 비교: padding 4

기본 GP는 100회, TV는 outer당 최대 200회·최대 5 outer이며,
수치 tolerance는 10⁻⁴, TV data discrepancy는 1%로 설정했다.
1%는 선언한 실용적 중단 기준이지 통계적으로 추정한 noise bound가 아니다.
이 padding에서는 정답 voxel 자체의 data mismatch가 1.67%이므로,
1% 이하로 맞추는 과정에 interpolation 오차의 보상도 포함될 수 있다.

| 방법 | RI 오차 | 관측 g 잔차 | 물리적 산란장 오차 | 내부 평균 Δn 편향 | 축방향 span |
|---|---:|---:|---:|---:|---:|
| Direct | 46.9056% | 18.0841% | 18.1251% | −10.4075% | 4.3 μm |
| GP | 36.5246% | 1.1863% | 1.1860% | −5.8848% | 8.0 μm, 경계 제한 |
| Lim TV | 21.6945% | 0.7214% | 0.7249% | −5.1009% | 5.6 μm |

Direct는 반복 수렴 판정의 대상이 아니며, GP와 TV는 이 cap에서 strict
수렴을 충족하지 못했다. TV의 최종 stationarity는 1.050 × 10⁻⁴로,
10⁻⁴ 기준을 조금 넘지만 임의로 통과 처리하지 않았다.
TV RI 오차의 내부/배경 성분은 각각 17.5028%와 12.8180%이다.

[단면 비교](exact_reconstruction.png),
[중심선·수렴 이력](exact_profiles_convergence.png),
[주파수 support와 missing cone](exact_missing_cone.png).

Direct와 GP에서는 중심값의 과대추정과 축방향 진동이 나타나며,
GP의 8 μm span에는 분리된 경계 artifact도 포함된다.
따라서 이를 단순한 균일한 축방향 elongation이라고만 부르면 부정확하다.
Direct에는 native Fourier deposition 자체의 근사 오차도 포함된다.
TV는 내부의 균질성을 개선하지만 축방향 경계의 흐림과 배경 편향이 남는다.

### Matched linear control

아래 표는 동일한 A로 만든 y = Aχ_true를 역산한 결과이며,
물리적 forward의 정확성을 추가 검증하는 표가 아니다.

| 방법 | RI 오차 | Matched g 잔차 | Exact 산란장과의 오차 | 수치 수렴 |
|---|---:|---:|---:|---|
| Direct | 46.6289% | 18.0362% | 19.0033% | 직접 변환 |
| GP | 35.8446% | 0.1462% | 1.6952% | 미충족 |
| Lim TV | 19.8241% | 0.1207% | 1.7024% | 미충족 |

TV의 stationarity는 1.155 × 10⁻⁴, axial span은 5.6 μm이다.
모델 불일치를 제거해도 현재 복원 오차가 남는다는 사실을 확인했지만,
유한 반복이므로 이 숫자를 missing cone의 불가피한 최소 오차라고
해석해서는 안 된다.
[Matched 단면](linear_control_reconstruction.png)과
[Matched 중심선·이력](linear_control_profiles_convergence.png)을 함께 저장했다.

## 7. Padding refinement와 수치 민감도

### 동일한 iteration cap에서 padding 8

아래는 기본 비교와 같은 GP 100회, TV 200 × 5회, ρ = ρ₊ = 0.1을
사용한 결과다. TV의 물리적 data/TV 비율은 유지했지만, splitting
penalty의 상대적 크기와 유한 반복의 진행률은 달라질 수 있다.

| 방법 | RI 오차 | 관측 g 잔차 | 물리적 산란장 오차 | 내부 평균 Δn 편향 | 수치 수렴 |
|---|---:|---:|---:|---:|---|
| Direct | 46.9056% | 17.5935% | 17.6344% | −10.4075% | 직접 변환 |
| GP | 36.3767% | 0.2221% | 0.2263% | −8.6089% | 100회 cap |
| Lim TV | 28.0123% | 0.3997% | 0.4134% | −10.6066% | 1000회/5 outer cap |

Direct는 native 격자 gridding을 사용하므로 padding을 바꾸어도 복원값은
같고, 재평가하는 forward의 interpolation이 달라져 field 잔차만 변한다.
GP의 작은 field 잔차에도 큰 RI 오차와 경계 artifact가 남는다.
TV의 stationarity는 2.231 × 10⁻⁴로 수렴 기준을 충족하지 못했다.
[Padding 8 단면](padding8_reconstruction.png)과
[중심선·이력](padding8_profiles_convergence.png)에 이 차이를 표시했다.

### TV weight, inner 반복, splitting penalty

아래 padding 4 실험은 outer당 최대 1000회, 최대 5 outer를 허용했다.
`alpha_relative`만 바꾸면 정규화 목적함수가 달라지고, ρ만 바꾸면 같은
inner 목적함수를 푸는 수치 경로가 달라진다. 유한 tolerance에서 outer에
전달되는 해와 reinjection 데이터도 달라질 수 있으므로 결과가 반드시
일치하는 것은 아니다.

| Run | alpha_relative | ρ = ρ₊ | RI 오차 | 관측 g 잔차 | 산란장 오차 | Inner 합계 / outer |
|---|---:|---:|---:|---:|---:|---:|
| tv_half | 2 | 0.1 | 28.1528% | 0.8185% | 0.8314% | 690 / 1 |
| tv_long | 4 | 0.1 | 25.9338% | 0.8419% | 0.8564% | 864 / 2 |
| tv_double | 8 | 0.1 | 27.7812% | 0.9694% | 1.0015% | 1019 / 2 |
| tv_penalty | 4 | 0.01 | 19.5501% | 0.7395% | 0.7469% | 421 / 2 |

네 경우 모두 선언한 10⁻⁴ inner 기준과 1% data discrepancy를 충족했다.
다만 `tv_long`은 capped 기본 TV보다 RI 오차가 커졌고, 같은 목적함수의
`tv_penalty`와도 6.38 percentage points 차이가 났다. 이는 더 오래
반복하거나 수치 수렴 flag를 얻는 것만으로 RI 정확도와 수치 안정성이
확립되지 않음을 보여준다. 정답 오차가 가장 작은 ρ를 골라 검증된 성능으로
제시하지 않고, 이 차이를 민감도 결과로 보존했다.
`tv_double`의 중심선 span은 8 μm로 경계에 닿으므로 두께 복원 성공이 아니다.
그 산란장 오차가 1%를 조금 넘는 것도 모순이 아니다. 중단 조건은 물리적
산란장이 아니라 선형화한 Ewald 데이터의 잔차에 적용된다.

### 물리적으로 penalty를 맞춘 padding 8

연산자 정규화 c²는 padding 4에서 0.01262063435, padding 8에서
0.00399053009이다. 따라서 ρ₈ = 0.01 × c₈²/c₄² = 0.00316190928을
사용해 data, TV와 splitting penalty의 상대적 크기를 맞췄다.
별도 회귀 시험에서 A와 y를 √s배, α와 두 ρ를 s배 했을 때 같은
20회 ADMM iterate를 10⁻¹⁰ 이내로 재현했다. Padding 변경은 여기에
interpolation 변경도 포함하므로 이 불변성 시험과 같은 문제는 아니다.

`tv_refined`는 240 inner / 1 outer에서 수치 수렴했지만 RI 오차는
29.3580%, data 잔차는 0.8719%, 산란장 오차는 0.9310%였다.
첫 inner 해가 1% discrepancy를 만족했으므로 이 run에서는 Lim 절차의
추가 data reinjection이 실제로 실행되지 않았다.
내부 평균 Δn는 15.45% 낮고, 중심선 half-contrast 영역은 계산 상자의
축방향 경계에 닿았다. 같은 문제의 tolerance를 10⁻⁵로 낮춘
`tv_strict`까지 완료한 비교는 다음과 같다.

| Run | Tolerance | RI 오차 | 관측 g 잔차 | 산란장 오차 | Inner / outer | Stationarity |
|---|---:|---:|---:|---:|---:|---:|
| tv_refined | 10⁻⁴ | 29.3580% | 0.8719% | 0.9310% | 240 / 1 | 1.847 × 10⁻⁵ |
| tv_strict | 10⁻⁵ | 28.9450% | 0.8689% | 0.9293% | 856 / 1 | 1.718 × 10⁻⁶ |

두 경우 모두 수치 수렴했고, strict run의 PCG 실패는 0회였다.
Tolerance를 10배 엄격하게 했을 때 RI 오차 변화는 0.4130 percentage
points이며, 두 복원값 자체의 contrast-normalized 차이는 1.6580%였다.
따라서 이 padding 8 설정의 약 29% RI 오차를 단순히 기존 inner 반복
수가 부족했기 때문이라고 설명할 수는 없다. 다만 이 두 tolerance의
비교만으로 정확한 최적해나 RI grid convergence를 증명한 것은 아니다.

Strict run의 내부/배경 오차 성분은 각각 17.0871%, 23.3632%이고,
내부 평균 n은 1.3356326233, 평균 Δn 편향은 −15.2254%이다.
전체 voxel의 최대 Δn도 정답의 87.87%이며 중심선 span은 여전히
8 μm·경계 제한이다. 이는 작은 field 잔차와 별개로, 굴절률의 과소평가와
배경으로 퍼진 구조가 남았다는 뜻이다. 첫 outer에서 종료되어 추가
data reinjection은 없었다.
[Strict TV 단면·중심선](tv_strict_reconstruction.png)과
[전체 TV 민감도 비교](tv_sensitivity.png)를 함께 저장했다.

### GP projection 정확도 검사

GP의 outer cap을 300회로 늘리고 CGLS tolerance를 10⁻⁵,
inner cap을 200회로 변경했다. Matched control은 158회에
iterate-change/normal-residual 기준을 만족했지만 RI 오차는 33.0805%,
matched 잔차는 0.1569%였다. 중심선 span 6.1 μm에도 분리된 경계
artifact가 포함되므로 이 길이만으로 축방향 복원이 정확하다고 볼 수 없다.
독립 exact 데이터의 동일 설정인 `gp_tight`도 300회 cap까지 완료했다.

| Run | 데이터 | RI 오차 | 해당 데이터 잔차 | Exact 산란장 오차 | 반복 / 수치 수렴 |
|---|---|---:|---:|---:|---|
| gp_control_tight | Matched | 33.0805% | 0.1569% | 1.6877% | 158 / 충족 |
| gp_tight | Exact | 85.0047% | 4.6220% | 4.6115% | 300 / 미충족 |

`gp_tight`의 최종 상대 iterate change는 0.01335로 10⁻⁴ 기준보다
훨씬 크고, 300회 중 189회의 CGLS correction이 선언한 normal-residual
기준을 충족하지 못했다. 마지막 correction만 기준을 통과했다고 전체
GP가 수렴한 것은 아니다. Matched control과의 차이는 모델 불일치와
ill-conditioning에 대한 민감도를 보여주지만, 부정확한 내부 projection도
섞여 있으므로 이 값을 이상적인 GP 최적해의 오차라고 해석할 수 없다.
단순히 CGLS tolerance와 iteration cap을 더 엄격하게 설정하는 것은
이 물리적 데이터에서 신뢰할 만한 개선책이 아니었다.

## 8. 해석상 제한과 다음 수정 후보

1. **Missing cone과 prior:** 조명·검출 NA를 넓혀도 현재 각도 집합에는
   축방향 미관측 정보가 남는다. 굴절률 비음수성과 TV만으로 정답
   spheroid를 유일하게 결정했다는 증거는 없다.
2. **TV 경계 조건:** periodic 차분은 양 끝 voxel을 이웃으로 취급한다.
   축방향 배경이 확장되거나 wraparound 형태가 생기는 데 영향을 줄 수
   있으며, 현재 경계 artifact의 원인을 이것 하나로 분리 검증하지는 않았다.
   물체 바깥 배경이 알려진 경우의 정당한 background 조건, 더 큰 복원
   상자, 비주기적 TV 경계는 후속 비교 후보다.
3. **GP:** 불일치 데이터에서 강한 data projection과 positivity의 충돌,
   CGLS의 불완전한 projection을 구분해야 한다. 이 문제에는 비음수
   least-squares나 오차를 허용하는 data-set projection도 대안이지만,
   이번 GP 결과를 그런 다른 알고리즘의 결과라고 부르지는 않는다.
4. **측정/모델:** 실제 잡음과 모델 오차에 맞춘 data fidelity 및 추가
   관측각·시료 회전은 후속 대상이다. 균질 spheroid가 확실한 prior라면
   a, b, nₚ의 parametric fitting은 더 작은 inverse 문제지만, 현재의
   unrestricted 3D RI 복원과는 다른 문제다.

이 대안은 원인을 분리하기 위한 후속 실험 제안이며, 이번 결과에 정답
support를 뒤늦게 넣거나 추가 알고리즘을 섞어 성능을 높이지 않았다.

## 9. 검증 및 재현

현재 코드·명령과 reference는 [README](../README.md), 시험 이력과 별도
forward suite의 제한은 [verification](../verification.md)에 있다.
모든 결과의 settings, iteration history, convergence flag를 MAT/CSV에
저장하며, Python으로 RI 오차를 독립적으로 재계산한다.
[Source/data SHA-256 목록](source_data.sha256)은 최종 inverse 코드,
재사용한 forward helper, 선택한 exact reference와 관측 cache를 식별한다.
[최종 artifact 검사](artifact_check.log)에서 11개 run·17개 복원의
설정·이력·수렴 flag와 11개 그림을 확인했다. Direct 3개를 제외한
14개 반복 복원 중 7개는 선언한 수치 기준을 충족했고, 7개는 미충족이다.
[최종 inverse 시험](verification_final.log)은 통과했고, MATLAB Code
Analyzer는 12개 MATLAB 파일에서 0건을 보고했다. 별도 forward suite는
`test_born_rytov.m:178`의 tolerance 유효성 검사 불일치에서 실패했으며,
그 기존 함수/시험은 수정하지 않아 repository 전체 통과를 주장하지 않는다.
`provisional/`은 교체된 초기 구현과 중단된 시험을 보존한 곳으로 최종
결과에 합치지 않는다. 코드 commit/push와 다른 section 수정은 하지 않았다.
실행시간에는 CPU 동시 공유, thread 수, 수학적으로 동등한 FFT 최적화의
차이가 있으므로 정밀한 알고리즘 속도 benchmark로 사용하지 않는다.

이번 데이터에는 잡음을 추가하지 않았다. Δn를 줄이면 모델 오차는
줄어들지만 산란 신호도 약해지므로, 본 결과는 고정 실험 잡음에서의
SNR 또는 실제 장비의 복원 성공률을 보장하지 않는다.
