# forward 개발·실험·오차 분석 이력

작성 기준: 2026-09-12  
범위: 이 대화에서 다룬 Born/Rytov forward 모델, 기존 prolate 결과, oblate 구조 재평가, xz 전기장 그림, Q의 정의·분포·물리적 해석까지.

## 1. 문서의 목적과 증거의 구분

이 문서는 결과 숫자뿐 아니라 질문이 어떻게 구체화되었는지, 어떤 정의를 유지하거나 수정했는지, 어떤 계산을 실행했고 어디에 저장했는지, 무엇을 검증했고 무엇을 아직 단정할 수 없는지를 보존한다. 수식은 사용자의 요청에 따라 LaTeX 명령 대신 일반 문자와 수학 기호로 쓴다.

기록은 다음 세 종류의 증거를 구분한다.

1. **이전 단계의 보존 기록:** 기존 [forward 설계](plan.md), [prolate detector 분석](tests/rytov_analysis.md), [prolate 원인 분석](exp/report.md), [spheroid solver 이력](../spheroid-analytic-forward/history.md)에 남아 있던 구현·실험.
2. **이번 oblate 평가의 실행 증거:** 재현 가능한 MATLAB script, CSV, MAT, 완료 marker를 포함한 log.
3. **후속 해석:** 저장된 exact 장에서 Q와 공간 주파수를 분석하고 광선 추적을 대조한 결과. 관측된 수치와 물리적 추론을 구분한다.

문서 작성 과정에서는 기존 자료와 현재 코드·결과 파일을 읽고 수치 및 경로를 대조했다. 과거의 모든 무거운 MATLAB 계산을 이 문서 작성 때문에 다시 실행한 것은 아니다. 과거 기록에만 있는 console 수치는 그 사실을 별도로 표시한다.

여기서 exact는 무한 급수를 기호적으로 닫힌 형태로 계산했다는 뜻이 아니라, spheroidal expansion을 충분히 수렴시키고 경계조건·cutoff·공개 evaluator 등과 대조한 수치적 full-wave reference를 뜻한다. Scalar exact와 Maxwell exact는 서로 다른 방정식의 해다.

## 2. 질문과 작업이 진행된 순서

| 단계 | 요청·쟁점 | 수행한 작업 또는 정리 |
|---|---|---|
| 1 | spheroid-analytic-forward의 exact 해와 과거 forward 구조의 관계 | 기존 reference와 forward 모델의 물리 조건을 확인했다. |
| 2 | prolate의 큰 Rytov 오차를 줄이기 위한 구조 변경 | 굴절률 대비를 낮추고 z 방향으로 납작한 oblate를 새 평가 대상으로 삼았다. |
| 3 | illumination NA=0.5, detector NA=1.0에서 scalar approximation 가능 여부 | 이전 scalar/Maxwell 비교의 의미와 허용 오차를 확인했다. |
| 4 | 회전대칭이면 방위각을 생략해도 되는가 | 원래 x polarization의 회전 convention을 직접 확인했다. |
| 5 | 같은 편광 정의로 Born/Rytov를 다시 평가 | 81개 입사 방향과 수치 refinement를 포함한 oblate 실험을 작성·실행했다. |
| 6 | Rytov가 약 10% 틀리는 원인 | scalar reduction, 근사 절단, 수치 오차를 구분하고 2차 로그항 및 대비 대조 계산을 수행했다. |
| 7 | exact/Born/Rytov의 공간 그림 | 최초의 3차원 그림 요청은 이후 y=0인 xz 단면의 전체장 크기 그림으로 구체화되었다. |
| 8 | 0°, 5°, 10°, 15°의 4×3 그림 | 행은 입사각, 열은 Exact/Born/Rytov인 공통 색상척도 그림을 저장했다. |
| 9 | 왜 10%가 5%가 되었는가 | 복소 산란장·복소 전체장·전체장 크기의 서로 다른 오차 정의를 구분했다. |
| 10 | Q의 정확한 식과 파장과의 관계 | 복소켤레 없는 내적, 진폭·위상 분해, 단위, Q와 f의 비교를 설명했다. |
| 11 | Q가 큰 곳이 꼭지점인가 | exact 배열에서 최대 위치와 분포를 계산하고 Q colormap을 추가했다. |
| 12 | 왜 내부 중앙에서 Q가 큰가 | 기울기 항 분해, 음의 kz 성분 확인, Snell 굴절과 한 번의 내부 반사 광선 추적으로 해석했다. |
| 13 | 모든 내용을 자세히 보존 | 이 파일 forward/history.md를 새로 작성했다. |

최종적으로 생성한 공간 그림은 굴절률 복원 결과나 3D isosurface가 아니다. 사용자가 마지막에 지정한 전체 Ex의 xz 단면이다. 서로 다른 입사각을 동시에 합친 focused beam도 구현하지 않았다.

## 3. 공통 물리 기호와 convention

| 기호 | 의미 |
|---|---|
| λ₀ | 진공 파장 |
| nₘ, nₚ | 배경 및 물체 굴절률 |
| k₀ = 2π/λ₀ | 진공 파수 |
| kₘ = nₘk₀, kₚ = nₚk₀ | 배경 및 물체 내부 파수 |
| a | z 방향 반축 |
| b | x, y 방향 반축 |
| θ, φ | 입사 방향의 polar angle과 azimuth |
| U₀ | scalar 입사장 |
| U | scalar 전체장 |
| Uₛ = U − U₀ | scalar 산란장 |
| Uᴮₛ, Uᴿₛ | Born 및 Rytov 산란장 |
| E, Eₛ | Maxwell 전체장 및 산란장 벡터 |
| e₀ | 단위 입사 편광 벡터 |
| P | detector의 이상적인 원형 pupil operator |
| ψ | 입사장으로 나눈 전체장의 복소 로그 |
| Q | ∇ψ · ∇ψ, 1차 Rytov에서 생략하는 항 |

시간 convention은 exp(−iωt), outgoing Green 함수는 G(R)=exp(ikₘR)/(4πR)이다. 중심이 원점인 spheroid는 다음과 같다.

(x/b)² + (y/b)² + (z/a)² < 1

산란 퍼텐셜은 다음과 같다.

f(r) = k₀²[n(r)² − nₘ²]

균질 물체 내부에서는 상수 f₀=k₀²(nₚ²−nₘ²), 외부에서는 f=0이다. 모든 길이는 같은 단위를 사용하며 이번 수치표에서는 μm를 사용한다.

Scalar Helmholtz 방정식은 다음과 같다.

(∇² + kₘ² + f)U = 0

재료는 이번 평가에서 균질·등방성·비흡수성이다. Q의 비선형성은 선형 파동방정식을 로그장으로 바꾸면서 생기는 비선형성으로, Kerr 효과 같은 비선형 광학 물질 응답을 넣은 것이 아니다.

## 4. 기존 forward 구현

### 4.1 연속 spheroid의 표면 적분 모델

[born_rytov_spheroid.m](born_rytov_spheroid.m)은 균질 prolate·oblate·sphere에 대해 임의의 내부·외부 점에서 scalar Born/Rytov 산란장을 계산한다.

Born 산란장은 다음 volume potential이다.

Uᴮₛ(r) = ∫ᵥ G(r−r′) f(r′) U₀(r′) dV′

균질 spheroid에서는 Green의 두 번째 항등식과 보조함수 w(r′)=(k̂ᵢ·r′)U₀(r′)/(2ikₘ)를 이용하여 이를 surface integral과 내부 jump term으로 바꾼다.

Uᴮₛ(r) = f₀{∮ₛ[G ∂ₙw − w ∂ₙG]dS′ − χᵥ(r)w(r)}

χᵥ는 내부에서 1, 외부에서 0이다. 따라서 내부 관측점의 volume-kernel singularity를 직접 3차원 적분하지 않는다. 표면에 정확히 놓인 점은 별도의 principal-value 처리가 필요하므로 현재 공개 함수는 offset을 요구한다.

주요 구현 결정은 다음과 같다.

- 표면 θ 적분은 Gauss–Legendre, φ 적분은 주기적 trapezoid이며 Nφ=2Nθ다.
- coarse/fine quadrature 쌍에서 Born과 Rytov를 각각 검사한다.
- Rytov에는 exp(ψ)−1 대신 복소 expm1을 사용한다.
- 멀리 떨어진 점에서는 공통 exp(ikₘr)을 분리하고 R−r를 cancellation-free 식으로 계산한다.
- 기본 tol=1e−8, Nmax=1200이며, 마지막 cap pair까지 실행한다.
- 표면 근처에서는 전역 quadrature refinement의 비용이 커질 수 있다.
- 수렴 실패나 nonfinite 출력은 기본적으로 error이며, warn 모드에서는 실패 mask와 validated=false를 남긴다.
- 이 함수의 수치 적분이 수렴했다는 사실은 Born/Rytov 근사 자체가 Maxwell 해와 가깝다는 뜻이 아니다.

공개 scalar forward 함수에는 polarization 인자가 없다. 입사 방향만 scalar 평면파에 들어가며, vector reference와 비교할 때 e₀ 방향 projection을 별도로 정의한다.

### 4.2 Voxel FFT/Ewald 모델

[born_rytov_voxel.m](born_rytov_voxel.m)은 3차원 굴절률 배열 n(ix,iy,iz)을 입력받아 +z 쪽 detector 평면의 Born/Rytov 산란장을 계산한다.

주요 순서는 다음과 같다.

1. f=k₀²(n²−nₘ²)를 만들고 zero padding한다.
2. 3D FFT에 dx³를 곱해 연속 Fourier transform convention에 맞춘다.
3. 각 입사 방향에 대해 q=kₛ−kᵢ인 shifted Ewald hemisphere에서 물체 spectrum을 보간한다.
4. kz>0인 propagating hemisphere에 i exp(ikz zdet)/(2kz)를 곱한다.
5. 2D inverse FFT로 detector 전의 linear Born field를 만든다.
6. 아래 절의 올바른 Rytov/pupil 순서로 최종 산란장을 만든다.

현재 보간은 linear interpn이며 기본 padding factor는 2다. Evanescent 및 정확히 grazing인 mode는 제외한다. Detector와 voxel volume은 같은 횡방향 pitch/FOV를 사용한다.

입력 검사는 detector Nyquist 조건뿐 아니라 pupil 적용 전 전체 propagating hemisphere와 입사 이동된 Ewald sample의 FFT 범위를 확인한다. 입사 reference가 detector pupil 안에 있고 zdet가 voxel grid의 +z 바깥에 있어야 한다. 따라서 작은 detector NA만 보고 voxel pitch를 크게 잡을 수 없다.

### 4.3 이전 Rytov detector-NA 순서 오류와 수정

이전 구현은 detector NA로 잘라낸 Born 장으로 Rytov 로그항을 만들었고, 지수화로 다시 생성된 pupil 밖 주파수도 최종 출력에 남겼다. 선형 pupil과 비선형 지수 연산은 교환되지 않는다.

수정된 순서는 다음과 같다.

Bpre = pupil을 적용하지 않은 forward-propagating Born 산란장  
ψ₁ = Bpre/U₀  
Uᴿₛ = P[U₀ expm1(ψ₁)]  
Uᴮₛ = P[Bpre]

따라서 info.rytov_phase는 pupil 이전의 ψ₁이며, 일반적으로 pupil 이후 log(1+Uᴿₛ/U₀)와 같지 않다. Forward 지수화에는 phase unwrapping이 필요하지 않지만, 측정장에서 로그를 역으로 구할 때의 branch/unwrap 문제는 별개다.

수정 후 sampling guard에는 dx<λ₀/(2nₘ)가 추가되었고, 이번 파장과 배경에서는 상한이 약 0.199194 μm다. dx=0.1 μm는 이를 만족한다. 입사각별 FFT bandwidth 검사도 따로 수행한다.

기존 보존 log의 NA regression은 Born 약 2.194e−16, Rytov 약 2.112e−16, pupil 밖 잔여 약 1.856e−16을 기록한다. 이는 operator 구현 검사이지 물리 근사 오차가 machine precision이라는 뜻이 아니다. 근거는 [이전 NA 수정 기록](tests/rytov_analysis.md)과 [verification.log](exp/verification.log)다.

## 5. 이전 prolate 결과: 새 oblate 실험의 출발점

### 5.1 Baseline과 중요한 수치

| 항목 | 기존 prolate |
|---|---:|
| λ₀ | 0.532 μm |
| nₘ / nₚ | 1.335381534 / 1.37 |
| a / b | 5 / 2.5 μm |
| 전체 x, y, z 길이 | 5, 5, 10 μm |
| 입사 / 편광 | +z / x |
| detector z / NA | 7 μm / 0.1 |
| voxel baseline | 128³, dx=0.1 μm, padding 2 |
| baseline FOV | 12.8×12.8 μm |
| 중심 ray 위상 지연 | 4.088613 rad |

같은 baseline의 complex scattered-field L2 error는 Born 197.560%, 수정된 voxel Rytov 25.735%였다. 수정 전 Rytov를 그대로 비교하면 74.396%, 지수화 후 pupil만 다시 씌운 값은 26.214%였으므로, 단순 재필터링과 올바른 operator 수정도 구분해야 한다. 더 이전 대화의 약 187%/33%는 단일 spectral line의 수치로, 전체 circular pupil의 결과와 같지 않다.

보존된 15-row 분석에는 padding 4, dx=0.05, FOV 확대, 633 nm, 물체 크기 절반, detector 거리, 낮은 굴절률 대비 control이 있다. 전체 원자료는 [rytov_analysis.csv](tests/rytov_analysis.csv)에 있다.

### 5.2 Scalar reduction과 Rytov 절단의 분리

별도의 experimental scalar prolate boundary-value reference는 같은 특수함수들을 사용하되 U와 ∂ₙU의 연속성을 만족시킨다. 이 코드는 real-index prolate의 축방향 scalar 입사에 제한되며, 후속 oblate batch가 추가되었다고 이 예전 API의 지원 범위가 자동으로 넓어진 것은 아니다.

- Scalar exact와 Maxwell Ex의 차이: 0.252653%.
- z=7 μm에서 연속형상 Rytov와 scalar exact의 차이: 24.1742%.
- 같은 연속형상 Rytov와 Maxwell Ex의 차이: 24.3021%.

따라서 이 baseline의 큰 잔차는 주로 scalar polarization reduction 밖의 문제였다. 25.735% voxel 값과 24.1742% 연속형상 scalar-reference 값은 수치 구현과 reference가 다르므로 동일 값으로 취급하지 않는다.

### 5.3 이전 prolate의 원인 분석에서 확인한 것

[기존 원인 보고서](exp/report.md)는 다음 결과를 보존한다.

| Detector z (μm) | 1차 Rytov | 2차 로그항 diagnostic |
|---:|---:|---:|
| 5.2 | 16.8730% | 4.3934% |
| 7 | 24.1742% | 13.3812% |
| 14 | 47.2136% | 2704.996% |

2차 Rytov는 모든 조건에서 개선되는 일반적 해결책이 아니다. 특히 z=14 μm의 실패를 새 oblate의 개선 결과로 지워서는 안 된다.

같은 중심 위상 지연 4.088613 rad와 detector 조건을 유지한 geometry control은 prolate 24.1742%, sphere 12.7176%, 무한 평판 2.9774%의 산란장 오차를 주었다. 평판의 전체장 오차 5.3004%는 분모가 다른 수치다. 중심 위상 지연 하나만으로 Rytov 정확도를 결정할 수 없다는 증거다.

또한 prolate exact scalar의 pupil 이전 전체장에서 ρ=0.1598748398 μm, z=6.996322478 μm의 ring zero를 찾았다. |U|는 L=112에서 2.41e−14, L=132에서 8.64e−12였다. 이 근처에서는 log(U/U₀)가 매끄럽고 유한하다고 가정할 수 없다. 이것은 NA-filtered detector에도 동일한 zero가 있다는 주장이나 이 한 ring이 전체 24%를 설명한다는 주장이 아니다.

예전 Q 분석은 3차원 체적 가중 regional RMS를 사용했고, 뒤의 oblate Q 표는 2차원 xz 단면 통계다. 두 종류의 RMS를 직접 비교하지 않는다.

## 6. 새 oblate 조건과 문헌·구조의 출처

사용자는 기존 index와 spheroid 설정의 배경으로 [ref/tomocube-rytov](../ref/tomocube-rytov)의 논문들을 언급했다. 해당 폴더에는 [Optics Express PDF](../ref/tomocube-rytov/oe-23-13-16933.pdf)와 [Scientific Reports PDF](../ref/tomocube-rytov/s41598-018-25886-8.pdf)가 있다.

다만 기존 보존 기록은 평균 핵 부피 약 300 μm³를 사용자가 제안한 배경값으로 구분하고, a=3 μm, b=5 μm를 논문의 특정 표에서 직접 추출한 축 길이라고 인증하지 않았다. 이 이력도 논문이 현재의 모든 숫자를 그대로 제시했다고 주장하지 않는다. 현재 형상은 문헌 배경을 참고하여 선택한 실험용 균질 모델이다.

| 항목 | 기존 prolate | 새 oblate |
|---|---:|---:|
| λ₀ | 0.532 μm | 동일 |
| nₘ | 1.335381534 | 동일 |
| nₚ | 1.370 | 1.365 |
| a, z 반축 | 5 μm | 3 μm |
| b, xy 반축 | 2.5 μm | 5 μm |
| 전체 x×y×z | 5×5×10 μm | 10×10×6 μm |
| detector z | 7 μm | 5 μm |
| +z 끝과 detector의 거리 | 2 μm | 2 μm |
| 주 비교 detector NA | 0.1 | 1.0 |
| illumination NA | 정상입사 중심 | 0–0.5 |
| 중심 ray 위상 지연 | 4.088613 rad | 2.09885 rad |

새 조건의 파생량은 다음과 같다.

- nₚ−nₘ = 0.029618466, 상대 굴절률 차이는 약 2.22%.
- λₘ=λ₀/nₘ ≈ 0.398388 μm.
- z 두께 6 μm는 약 15개의 배경 매질 파장에 해당한다.
- 부피 V=4πab²/3=314.159265 μm³.
- semifocal length √(b²−a²)=4 μm, oblate 표면 좌표 ξ₀=a/4=0.75.
- f₀=11.156402219375025 μm⁻².
- illumination NA=nₘ sinθ, 따라서 NA=0.5는 θ≈21.98877°다.

이 모델은 핵·세포질의 이중 경계, membrane, organelle, 흡수, 굴절률 분산을 포함하지 않는다. nₚ가 일정한 하나의 spheroid와 배경만 있다.

## 7. x polarization과 방위각에 대한 질문

### 7.1 원래 코드의 정의

[incident_plane_wave.m](../spheroid-analytic-forward/incident_plane_wave.m)의 pol=[1,0]을 그대로 유지했다.

k̂ᵢ = (sinθ cosφ, sinθ sinφ, cosθ)

eTE = (−sinφ, cosφ, 0)

eTM = (cosθ cosφ, cosθ sinφ, −sinθ)

e₀ = cosφ eTM − sinφ eTE = Rz(φ) Ry(θ) Rz(−φ) x̂

이는 정상입사의 x 방향을 입사 방향으로 추가 roll 없이 기울이는 convention이다. 방위각을 바꿀 때 언제나 eTM만 쓰는 co-rotating polarization과는 다르다.

특수한 경우는 다음과 같다.

- θ=0: e₀=x̂.
- φ=0: e₀=(cosθ,0,−sinθ).
- φ=90°: e₀=x̂.
- 모든 경우: e₀·k̂ᵢ=0.

입사각을 바꾸면 transverse 조건을 만족하도록 전기장 방향도 바뀐다. 반면 정상입사를 유지한 채 detector NA만 늘리는 실험에서는 입사 편광이 회전하는 것이 아니라 더 넓은 산란 방향을 수집한다.

### 7.2 회전대칭이면 φ를 생략할 수 있는가

축대칭 물체와 원형 pupil만으로 모든 polarization convention의 방위각 독립성을 자동으로 주장할 수는 없다. 원래 pol=[1,0]은 φ에 따라 TE/TM 혼합 비율이 변한다.

따라서 각 θ에서 TM과 TE를 모두 풀고, φ에 따른 회전과 cosφ, −sinφ 혼합으로 원래 x polarization을 복원했다. Voxel/연속형상 forward는 각 입사 방향을 같은 Cartesian FFT grid에서 직접 평가했다. 이는 grid와 interpolation의 방향 의존성도 검사하기 위한 것이다.

이번 최종 co-polarized Rytov metric에서는 같은 illumination NA에서 azimuth에 따른 최대 차이가 약 0.0091 percentage points로 작았다. 따라서 이번 조건에서는 영향이 작다고 말할 수 있지만, 처음부터 물체의 축대칭만으로 정확히 같다고 가정한 것은 아니다.

## 8. 앞선 oblate scalar/Maxwell pilot과 현재 실험의 차이

기존 [spheroid history의 13절](../spheroid-analytic-forward/history.md)은 먼저 Born/Rytov가 아닌 scalar exact와 Maxwell exact를 비교했다. 그 단계는 MATLAB session의 console 결과를 문서로 남긴 것으로, 당시의 독립 script/CSV/MAT는 저장되지 않았다고 명시되어 있다.

정상입사 pilot의 주요 결과는 다음과 같으며 모두 %다.

| Detector NA | 편광 외 성분 ε⊥ | Scalar co 오차 | Scalar full-vector 오차 |
|---:|---:|---:|---:|
| 0.10 | 1.902210 | 0.075638 | 1.903713 |
| 0.30 | 2.478660 | 0.174686 | 2.484804 |
| 0.50 | 2.715368 | 0.303327 | 2.732245 |
| 1.00 | 3.023760 | 0.726344 | 3.109697 |
| 1.20 | 3.154261 | 1.079322 | 3.333637 |
| 1.30 | 3.263472 | 1.468017 | 3.578132 |

이 pilot에서 full-vector 오차의 첫 1% crossing은 detector NA 약 0.0281, co 오차만의 첫 1% crossing은 약 1.1755였다. 즉 scalar approximation이 괜찮다는 말에는 어떤 field component와 허용 오차를 말하는지 반드시 따라야 한다.

이후 illumination NA=0:0.05:0.5와 추가 θ=15°, 총 12개 polar condition을 조사했고, φ=0°,45°,90° 및 TE/TM 원형 pupil 적분의 대칭성을 확인했다. Detector NA=0.5–1.3의 801개 표본에서 조사한 최악의 scalar full-vector 오차는 최대 약 3.697%로 5% 안이었다. 정상입사 결과를 그대로 oblique의 작은 detector NA에 적용하면 안 된다. 고정된 +z pupil이 기울어진 강한 forward lobe를 놓치면 상대 오차가 커질 수 있다.

당시 시작했던 별도의 θ=15° 및 최대 입사각 public solve는 긴 validation 때문에 중단되었으므로, 그 두 실행의 validated=true를 얻었다고 기록하지 않는다. 대신 별도 batch 경계·cutoff·far-limit 검사로 결과를 확인했다. 이전의 illumination NA=0.9 언급은 사용자 정정에 따라 최종 0.5로 대체되었다.

현재 [oblate_na_sweep.m](exp/oblate_na_sweep.m)은 이 접근을 재현 가능한 실험으로 만들고, 실제 Born/Rytov 근사 오차와 수치 refinement까지 평가한 후속 작업이다. 앞선 console-only 12-condition pilot과 현재 81-direction sweep를 같은 실행으로 합치지 않는다.

## 9. 재현 가능한 oblate NA sweep의 구현

### 9.1 입사각 표본

최종 sweep는 illumination NA=0:0.05:0.5의 11개 반경과 φ=0:45:315°의 8개 방위각을 사용했다. 정상입사는 중복 계산하지 않아 총 1+10×8=81개 방향이다.

각 θ에서 incidence plane을 xz로 잡고 TM/TE 두 개의 independent right-hand side를 풀었다. Signed azimuthal order m마다 같은 Galerkin matrix를 공유하고 여러 θ와 편광을 동시에 풀어, 각 입사각마다 전체 solver를 반복하는 비용을 줄였다.

Batch reference는 다음 public 수학 구성요소를 재사용한다.

- sph_eigen: oblate angular eigenbasis.
- sph_angular: metric-weighted angular 함수.
- sph_radial: 내부 regular radial 함수.
- sph_radial_ode와 sph_radial_ode_eval: 외부 outgoing radial basis.
- sph_vecwave: vector M/N wave.
- gauss_legendre, incident_plane_wave: quadrature와 입사 convention.

Scalar system은 U와 ∂ₙU를, vector system은 Maxwell 경계조건을 사용한다. Column scaling 후 경계 matching을 풀며, solve grid와 다른 θ grid에서 residual을 평가한다. 재구성된 전체 경계장과 버린 incident azimuthal tail을 검사하므로, 거의 0인 개별 high-m RHS로 나눈 불안정한 sector 상대값을 최종 residual로 쓰지 않는다.

### 9.2 Cutoff와 reference 저장

Coarse (L,M)=(109,55), fine (114,59)를 사용했다. 입사장의 azimuthal quadrature는 256점이다. 각 batch의 여러 입사각과 TM/TE를 함께 평가한다.

[oblate_reference_0.mat](exp/oblate_reference_0.mat), [oblate_reference_1.mat](exp/oblate_reference_1.mat)은 far-field coefficient와 검증 정보를 보존한다. 저장 전 거대한 outgoing ODE dense trajectory인 rod_ext.solution을 제거했으므로, 이 compact MAT를 그대로 public near-field evaluator에 넣을 수는 없다.

원래 source comment에는 reference stage를 다시 실행해 full structure를 얻는다는 취지의 문구가 있지만, 실제 stage는 기존 compact checkpoint가 있으면 이를 읽는다. Near-field를 다시 계산할 때는 full batch를 새로 만드는 xz 또는 xz_validate stage를 사용해야 한다. Compact MAT의 dense trajectory가 자동 복원된다고 가정하지 않는다.

### 9.3 연속형상과 voxel의 분리

수치적 경계 voxelization과 Ewald interpolation 오차를 분리하기 위해 spheroid의 정확한 Fourier shape transform을 사용하는 control을 만들었다.

χ = √[b²(qx²+qy²)+a²qz²]

F(q) = f₀V × 3(sinχ−χ cosχ)/χ³

χ→0에서는 1−χ²/10+χ⁴/280의 series를 사용한다. 코드 일부가 이 shape argument를 Q라는 변수로 부르지만, 이후 로그장 비선형항 Q=∇ψ·∇ψ와는 다른 양이다. 이 문서에서는 shape argument를 χ로 구분한다.

이 연속형상 계산도 1차 Born/Rytov 방정식 자체는 그대로다. 정확한 형상 transform을 사용했다고 Rytov 근사가 exact가 되는 것은 아니다.

### 9.4 실행 configuration

| Configuration | 횡방향 N | dx (μm) | Nz | Padding | 입사 방향 수 |
|---|---:|---:|---:|---:|---:|
| baseline | 256 | 0.1 | 80 | 2 | 81 |
| window_51 | 512 | 0.1 | 연속형상 | 없음 | 81 |
| sampling_0p05 | 1024 | 0.05 | 연속형상 | 없음 | 81 |
| padding_3 | 256 | 0.1 | 80 | 3 | 7 |
| voxel_0p0667 | 384 | 1/15 | 120 | 2 | 7 |
| voxel_window_51 | 512 | 0.1 | 80 | 2 | 7 |

대표 7개 방향은 정상입사 1개와 illumination NA=0.25,0.5에서 φ=0°,45°,90°다. 모든 sweep row를 합친 [oblate_angles.csv](exp/oblate_angles.csv)는 264행이고, [oblate_summary.csv](exp/oblate_summary.csv)는 33개 configuration/metric summary다.

## 10. Detector 오차의 정의와 입사각 평균

### 10.1 Reference의 세 가지 비교

모든 detector field에 같은 이상적인 circular pupil을 적용한 뒤 비교한다.

u_co = Eₛ · e₀*

ε_co = ‖Umodel,ₛ − u_co‖₂ / ‖u_co‖₂

ε_scalar = ‖Umodel,ₛ − Uexact,ₛ‖₂ / ‖Uexact,ₛ‖₂

ε_full = ‖Umodel,ₛ e₀ − Eexact,ₛ‖₂ / ‖Eexact,ₛ‖₂

여기서 *는 complex conjugation이다. Polarization projection의 conjugation은 뒤의 Q에 없는 conjugation과 혼동하지 않는다. Scalar result에 e₀를 곱한 full-vector 비교는 scalar 모델의 vector reference 오차이며, 별도의 vector Born/vector Rytov 구현이 아니다.

Scalar approximation 자체의 오차도 Umodel 대신 scalar exact를 넣어 측정한다. Polarization leakage와 co-polarized amplitude/phase mismatch는 다른 성분이다.

원형 pupil의 반경은 k₀ NA_det다. NA=1.0은 수집 범위이며, sample의 각 점에 적용하는 공간 cutoff가 아니다.

### 10.2 Far amplitude에서 detector spectrum으로의 변환

Eₛ(r) ≈ A(k̂ₛ) exp(ikₘr)/r

S(kx,ky;zdet) = 2πi A(k̂ₛ) exp(ikz zdet)/kz

이 변환과 Parseval 관계로 같은 detector spectrum에서 complex-field L2 norm을 비교한다. Reference cutoff 검사의 radial NA 적분 weight는 공통 상수를 제외하면 NA/(nₘ²−NA²)다. 이는 아래 illumination 평균의 annulus-area weight와 다르다.

실제 objective의 polarization transport, aberration, analyzer, apodization은 넣지 않았다. 유한 FFT window와 실제 camera crop도 같은 뜻은 아니다.

### 10.3 Illumination 평균

입사각별 오차를 단순히 11개 polar ring에 똑같이 배분하지 않았다. Illumination pupil의 면적을 균일하게 덮는다는 가정을 사용했다.

인접 NA 표본의 중간점을 annulus 경계로 하고 [0,0.5]에서 자른다. 각 annulus의 면적 비율을 같은 ring의 azimuth 표본에 균등 배분한다.

ε_mean = Σⱼ wⱼ εⱼ, Σⱼ wⱼ=1

ε_pooled = √[Σⱼ wⱼ ‖ΔUⱼ‖₂² / Σⱼ wⱼ ‖Uref,ⱼ‖₂²]

Weighted mean과 pooled L2는 서로 다른 평균이다. 서로 다른 입사장들을 coherent하게 더한 뒤 오차를 계산하지 않는다. 연속 polar-angle 전 범위의 엄밀한 bound가 아니라 명시된 표본과 가중치의 결과다.

## 11. 새 oblate의 Born/Rytov 성능

### 11.1 주 결과

아래는 co-polarized Maxwell reference에 대한 complex scattered-field 오차이며 모든 값은 %다. 연속형상 최종 값은 51.2 μm window, dx=0.05 μm의 refined configuration이다.

| 모델 | Pupil-area weighted mean | Pooled L2 | 최소 | 최대 |
|---|---:|---:|---:|---:|
| 기존 voxel Born, baseline | 89.679 | 89.716 | 88.095 | 91.968 |
| 기존 voxel Rytov, baseline | 9.751 | 9.775 | 8.500 | 10.844 |
| 연속형상 Born | 91.319 | 91.368 | 89.084 | 93.648 |
| 연속형상 Rytov | 10.509 | 10.528 | 9.869 | 11.192 |
| Scalar exact | 0.705 | 0.705 | 0.658 | 0.733 |

주 결론은 다음과 같다.

- Born은 현재의 낮아진 대비와 납작한 구조에서도 약 91%로 부정확하다.
- 1차 Rytov는 크게 개선되지만 약 10.5%를 남긴다.
- Scalar exact와 Maxwell co-polarized field의 차이는 약 0.70%다.
- 따라서 scalar approximation의 성립을 근거로 Rytov도 같은 수준으로 정확하다고 판단할 수 없다.

Scalar exact를 reference로 한 연속형상 평균은 Born 91.281%, Rytov 10.525%였다. Maxwell co reference의 91.319%,10.509%와 거의 같으므로, 이 조건의 주요 잔차는 scalar polarization reduction으로 설명되지 않는다.

Scalar field를 e₀ 방향 vector로 만들어 Maxwell 전체 산란장과 비교하면 연속형상 Born/Rytov 평균은 91.327%/10.935%, baseline voxel은 89.689%/10.210%다. Scalar exact 자체의 full-vector 평균은 3.117%, 최대는 약 3.166%다.

### 11.2 입사각 의존성

아래는 각 illumination NA의 azimuth 평균이며 %다.

| Illumination NA | 연속형상 Born co 오차 | 연속형상 Rytov co 오차 |
|---:|---:|---:|
| 0.0 | 89.084 | 9.869 |
| 0.1 | 89.252 | 9.921 |
| 0.2 | 89.766 | 10.065 |
| 0.3 | 90.646 | 10.310 |
| 0.4 | 91.923 | 10.680 |
| 0.5 | 93.648 | 11.187 |

Rytov의 최악 표본은 illumination NA=0.5, φ=0°의 약 11.192%다. 같은 NA에서 방위각 차이는 최종 연속형상 co metric에서 최대 약 0.0091 percentage points였다.

기존 prolate의 25.74%와 현재 oblate의 약 10.5%는 geometry, index, illumination, detector NA와 z가 함께 달라진 결과다. 이 차이를 NA를 늘린 효과 하나 또는 shape 효과 하나로 분리하지 않았다.

## 12. Numerical audit와 수치 오차의 한계

### 12.1 Reference 검증

| 검사 | 상대 차이 또는 residual |
|---|---:|
| Vector cutoff (109,55) 대 (114,59) | 6.5277e−13 |
| Scalar cutoff | 6.7525e−13 |
| 독립 grid의 최대 vector boundary residual | 4.6597e−12 |
| 독립 grid의 최대 scalar boundary residual | 1.9486e−12 |
| Analytic far limit 대 public evaluator의 외삽 | 4.1261e−8 |

Boundary residual은 surface-weighted L2이지 최대 pointwise residual이 아니다. Batch에는 별도의 검증이 있으며, public spheroid_solve의 모든 adaptive runtime-validation contract를 그대로 통과했다는 뜻으로 validated=true를 붙이지 않았다.

작은 크기·낮은 대비의 구현 검사에서는 scalar Born amplitude와 transverse-projected vector Born amplitude를 각각 비교했다. Original polarization의 φ=90° 조건도 assertion으로 확인했다.

### 12.2 연속형상 grid와 FOV

| 변화 | 관측한 최대 영향 |
|---|---:|
| Window 25.6→51.2 μm | Rytov co 오차값 변화 약 0.0049 percentage points 이하 |
| 같은 물리 주파수에서 25.6→51.2 μm field 비교 | 약 0.0739% |
| 같은 물리 주파수에서 51.2→102.4 μm field 비교 | 약 0.0337% |
| 51.2 μm 고정, dx=0.1→0.05 μm field 비교 | 약 0.000101% |
| Illumination NA step 0.05→0.1 | 평균 오차 차이 0.0110 percentage points |
| Azimuth step 45°→90° | 평균 오차 차이 0.000048 percentage points |

즉 현재 약 10.5%라는 연속형상 Rytov 오차를 단순 FFT sampling 부족으로 설명할 수는 없다. 다만 mode residual의 10⁻¹² 수준을 최종 forward observable의 절대 정확도 보증으로 해석하지 않는다.

### 12.3 Voxel interpolation과 오차 상쇄

Baseline padding 2에서는 voxel Rytov와 연속형상 Rytov의 field norm 차이가 대표 방향에서 최대 약 4.505%였다. Padding 3에서는 최대 약 2.027%로 줄었으며, voxel의 Maxwell-relative error 값 자체는 최대 0.748 percentage points 바뀌었다.

Voxel pitch를 0.1→1/15 μm로 줄이는 control의 오차값 변화는 최대 약 0.0363 percentage points, detector window를 51.2 μm로 늘리는 control은 약 0.0847 percentage points였다.

Baseline voxel 평균 9.751%가 연속형상 10.509%보다 작다고 voxel Rytov가 물리적으로 더 정확한 것은 아니다. Numerical error와 approximation error가 부분적으로 상쇄된다. Voxel 구현의 field 자체가 1% 이하로 수렴했다고 증명한 것은 아니다.

전체 수치표와 해석은 [oblate_report.md](exp/oblate_report.md)에 보존되어 있다.

## 13. 0°, 5°, 10°, 15°의 Ex colormap

### 13.1 최종 그림의 정의

사용자가 최종적으로 지정한 그림은 다음과 같다.

- y=0의 xz 평면.
- +z로부터 θ=0°,5°,10°,15° 입사, φ=0.
- 원래 pol=[1,0] 유지.
- 행 4개: 입사각.
- 열 3개: Exact, Born, Rytov.
- 색: 입사장을 포함한 전체 |Ex|/|E₀|.
- 모든 12개 panel에 동일한 색 범위 0–2.8.
- 흰 선: (x/5)²+(z/3)²=1인 실제 spheroid 단면.
- 화살표: 입사 진행 방향.
- x,z의 실제 길이 비율을 같은 aspect로 표시.

Exact는 Maxwell의 laboratory Ex이다. Scalar Born/Rytov는 다음 projection을 사용한다.

Ex,B = cosθ [U₀+Uᴮₛ]

Ex,R = cosθ U₀ exp(Uᴮₛ/U₀)

입사 전기장 vector의 크기 |E₀|=1을 유지했으며, 각 panel의 incident Ex를 별도로 1로 재정규화하지 않았다. θ=15°에서 incident Ex는 cos15°다.

### 13.2 관측 범위와 계산 방식

| 항목 | 값 |
|---|---:|
| x 범위 | −8–8 μm |
| z 범위 | −5–9 μm |
| Nx × Nz | 201 × 176 |
| dx, dz | 0.08 μm |
| 전체 평가점 | 35,376 |
| 입사각 | 0°, 5°, 10°, 15° |
| 해당 최대 illumination NA | 약 0.346 |
| reference cutoff | L=114, M=59 |

이 그림은 물체 내부와 바깥의 raw total field다. Detector NA=1 filter를 물체 내부 단면에 적용하지 않았다. Pupil을 사용한 detector 비교와 full Green-function 공간 그림을 구분한다.

35,376개 점에서 surface quadrature를 반복하는 대신, 산란 퍼텐셜의 t=0 근처 symmetric derivative로 full-Green Born field를 계산했다.

n(t)=√[nₘ²+t(nₚ²−nₘ²)]

Uᴮₛ ≈ [U(+h)−U(−h)]/(2h), h=0.005

외부에서는 공통 exterior basis의 coefficient를 먼저 차분할 수 있다. 내부에서는 n(t)에 따라 regular basis가 달라지므로 각각의 전체장을 평가한 후 차분했다. 외부의 incident term은 정확히 분리했고, 내부에서도 symmetric difference에서 공통 입사장이 상쇄된다.

이 방법은 h→0에서 first Born이 되는 성질을 이용한 수치 평가이며, 유한 h의 정확성은 기존 surface-integral Born 함수 및 h refinement로 별도 확인했다. Reference index를 1.365로 유지한 채 Maxwell 장을 Born으로 부른 것이 아니다.

### 13.3 저장과 렌더링

[oblate_xz_fields.mat](exp/oblate_xz_fields.mat)에 exact, bornEx, rytovEx, scalar, u0, born, Q, inside, x, z, theta_deg 및 cutoff 검증값을 저장했다. [oblate_xz_exact.mat](exp/oblate_xz_exact.mat)에는 Cartesian scalar gradient도 있다.

MATLAB의 meshgrid 결과를 (:)로 펼친 배열이므로 Python에서는 (Nz,Nx), order="F"로 복원했다. [plot_oblate_xz.py](exp/plot_oblate_xz.py)는 이 배열을 읽어 Matplotlib로 렌더링한다. 보간은 nearest이고, 공통 color limit은 모든 배열의 실제 최대값을 포함하도록 정했으며 percentile clipping이나 panel별 rescaling은 하지 않았다.

결과물:

- [4×3 PNG](exp/oblate_Ex_xz_4x3.png)
- [4×3 PDF](exp/oblate_Ex_xz_4x3.pdf)
- [xz 상세 보고서](exp/oblate_xz_report.md)

그림을 실제로 열어 열·행 순서, 축, 색상척도, boundary, 편광 projection, 배열 방향을 확인했다.

## 14. xz 계산 중 발견한 좌표 roundoff 오류

첫 near-field run은 x=y=0인 일부 점에서 sph_angular의 coordinate-range validation에 걸렸다.

오류는 sph_coords의 oblate 계산에서 eta_abs가 반올림 때문에 1보다 아주 조금 커진 것이었다. 원래 격자에서 z=4.04,7.16,7.24,7.8 μm인 네 축 점이 이 문제를 재현했다. 이는 물리적 singularity가 아니라 허용 범위를 넘는 floating-point coordinate 값이었다.

해결은 공유 coordinate 함수 [sph_coords.m](../spheroid-analytic-forward/sph_coords.m)에 다음 한 줄을 추가하는 것이었다.

eta_abs = min(1, eta_abs);

후속 angular function에서 결과를 숨기거나 각 caller에 workaround를 넣지 않고, 수학적 범위가 정해진 coordinate 생성 지점에서 수정했다. [test_sph_coords.m](../spheroid-analytic-forward/tests/test_sph_coords.m)에 동일한 axial grid의 |eta|≤1과 basis 검사를 추가했다.

직접 grid assertion으로 수정 전 실패를 확인했고, 수정 뒤 test_sph_coords와 test_coords_vecwave가 통과했다. 처음 실패한 xz run은 재시작하여 끝까지 완료했다. 실패 내용은 대화에 남아 있고 oblate_xz.log는 성공 재실행으로 덮였으므로, 현재 log 하나에 첫 실패까지 보존되어 있다고 주장하지 않는다.

초기 oblate detector-sweep 계획은 production code를 바꾸지 않는 실험이었다. 이후 raw xz plot 과정에서 위 shared coordinate fix가 추가되었으므로, 최종 전체 작업을 "어떤 production code도 수정하지 않았다"라고 요약하면 부정확하다. Born/Rytov production 수식은 이번 oblate 실험에서 바꾸지 않았다.

## 15. 10%와 5%가 혼동된 이유와 정정

### 15.1 서로 다른 세 가지 지표

사용자가 왜 10%에서 5%로 감소했는지 질문했다. 이는 같은 알고리즘·같은 측정량에서 정확도가 갑자기 개선된 것이 아니다.

| 비교 | 결과 | 의미 |
|---|---:|---|
| 기존 81-direction detector 비교 | 평균 10.509% | Pupil 이후 복소 산란장, Maxwell co reference |
| 새 xz 단면의 복소 전체 Ex | 8.275–9.000% | 같은 raw 공간 window에서 위상과 크기 모두 비교 |
| 새 xz 단면의 전체 ∣Ex∣ | 5.220–5.808% | 동일 xz 데이터의 크기만 비교 |

마지막 표 행의 ∣Ex∣는 절댓값 기호이며 intensity ∣Ex∣²가 아니다. 세부 정의는 다음과 같다.

ε_amp = ‖|Ex,model|−|Ex,exact|‖₂ / ‖|Ex,exact|‖₂

ε_complex = ‖Ex,model−Ex,exact‖₂ / ‖Ex,exact‖₂

같은 xz 배열에서 8–9%가 5–6%가 되는 것은 위상 차이를 평가에서 제거했기 때문이다. 하지만 10.509%와 xz 결과 사이에는 관측 위치·pupil·전체장/산란장 분모·입사각 표본 및 평균 방식도 달라지므로, 전체 차이를 위상 오차 하나로만 설명할 수 없다.

### 15.2 xz 오류표

모든 값은 표시된 전체 xz window에 대한 %다.

| θ | Born 크기 | Born 복소장 | Rytov 크기 | Rytov 복소장 | Scalar exact Ex 대 Maxwell Ex, 복소장 |
|---:|---:|---:|---:|---:|---:|
| 0° | 58.689 | 83.821 | 5.220 | 8.275 | 2.242 |
| 5° | 58.928 | 84.241 | 5.291 | 8.361 | 2.276 |
| 10° | 59.695 | 85.563 | 5.494 | 8.601 | 2.437 |
| 15° | 60.962 | 87.766 | 5.808 | 9.000 | 2.687 |

정상입사에서 window 내 최대 |Ex|는 Exact 1.470, Born 2.774, Rytov 1.590이었다. Born은 downstream amplitude를 크게 과대평가한다.

[oblate_xz_plot_errors.csv](exp/oblate_xz_plot_errors.csv)는 전체 window와 물체 내부를 따로 기록하여 총 24행이다. 위 표는 전체 window만 뽑은 것이다. Raw scalar exact/Maxwell 오차 2.24–2.69%를 기존 pupil 이후 co 산란장의 약 0.70%와 같은 metric으로 취급하지 않는다.

이 지표 변경을 첫 답변에서 더 분명히 구분했어야 한다는 점을 대화에서 정정했다.

## 16. Q의 정확한 정의와 사용자의 수식 표기 요청

### 16.1 정확한 로그장 방정식

U≠0인 영역에서 다음과 같이 정의한다.

ψ = ln(U/U₀)

U₀ = exp(ikᵢ·r)

∇²ψ + 2ikᵢ·∇ψ + Q = −f

Q = ∇ψ·∇ψ = (∂ψ/∂x)² + (∂ψ/∂y)² + (∂ψ/∂z)²

1차 Rytov는 Q를 생략한다. ∇²ψ는 남기므로 회절 자체를 버리는 근사는 아니다. Q는 복소켤레 없는 내적이므로 일반적으로 복소수이며, 양의 실수인 ∇ψ*·∇ψ와 다르다.

대화에서 사용자는 LaTeX 코드 대신 기호로 식을 써 달라고 요청했다. 이후 설명의 기본 표기는 위와 같은 일반 문자·Unicode 식으로 정리했다.

### 16.2 진폭과 위상으로 분해

U = A exp(iΦ), U₀ = A₀ exp(iΦ₀)

α = ln(A/A₀), ΔΦ = Φ−Φ₀, ψ = α+iΔΦ

Q = ∣∇α∣² − ∣∇ΔΦ∣² + 2i(∇α·∇ΔΦ)

여기서 α와 ΔΦ는 실수이므로 각 gradient의 norm square는 일반적인 실수 norm이다. Q에는 로그 진폭의 공간 변화, 상대 위상의 공간 변화, 두 변화의 교차항이 모두 포함된다. Q 자체의 크기가 작더라도 개별 기울기가 각각 작다는 뜻은 아니며 항 사이 상쇄가 가능하다.

실제 계산은 principal log를 직접 취하지 않고 exact scalar의 analytic Cartesian derivative를 사용했다.

∇ψ = (∇U)/U − ikᵢ

Q = [(∇U)/U − ikᵢ] · [(∇U)/U − ikᵢ]

따라서 위상의 2π branch나 unwrap 선택이 이 Q 계산에 직접 들어가지 않는다. U=0에서는 로그장과 이 식이 singular하므로 별도 주의가 필요하다.

### 16.3 파장보다 작아야 하는가

ψ는 무차원, ∇ψ의 단위는 길이⁻¹, Q와 f의 단위는 길이⁻²다. 따라서 Q를 파장 λ 자체와 직접 비교할 수 없다.

물체 내부에서 생략항을 평가하는 국소 비교는 다음과 같다.

∣Q∣ ≪ ∣f∣

f = (2π/λ₀)²(n²−nₘ²)

λ₀²∣Q∣ / [4π²∣n²−nₘ²∣] ≪ 1

현재 균질 물체에서는 ∣f∣≈11.1564 μm⁻²다. 따라서 단지 Q가 1보다 작다거나 λ보다 작다는 표현은 단위와 대비를 빠뜨린 것이다. λ를 바꾸면 exact 장과 Q도 함께 바뀌므로, 위 식만 보고 파장을 줄이면 무조건 좋아진다고 결론 내리지 않는다.

이 비율은 국소적인 생략항 진단이며 detector error 10%를 바로 계산하는 식이 아니다. Exterior에서는 f=0이므로 Q/f를 정의해 같은 방식으로 비교하지 않는다. 외부에서도 Q가 반드시 0인 것은 아니다.

### 16.4 현재 구조에서 작은가

현재 oblate의 exact scalar에서 Q가 모든 위치에서 충분히 작다는 가정은 성립하지 않았다. 많은 점에서는 작지만 일부 좁은 영역에서 커진다.

한편 중심 ray 위상 지연 2.1 rad 자체가 Rytov의 실패 조건은 아니다. 큰 누적 위상을 지수로 표현하는 것과, 로그장의 공간 기울기에 의한 고차 보정을 정확히 표현하는 것은 다른 문제다.

## 17. 2차 로그항과 대비 sweep으로 원인을 검사

### 17.1 독립 계수 추출

산란 퍼텐셜에 t를 곱하고, 형상·배경·파장을 유지했다.

n(t)=√[nₘ²+t(nₚ²−nₘ²)]

U(t)=U₀+tU₁+t²U₂+…

ψ₁=U₁/U₀

ψ₂=U₂/U₀−ψ₁²/2

Uᴿ¹(t)=U₀ exp(tψ₁)

Uᴿ²(t)=U₀ exp(tψ₁+t²ψ₂)

U₁과 U₂는 약한 양·음의 potential scale에서 얻은 exact scalar solution의 symmetric difference로 추출했다. 물리적 대비 t=1의 오차를 최소화하도록 fit한 보정항이 아니다.

Total field를 사용하면 U₂≈[U(+h)+U(−h)−2U₀]/(2h²)이고, scattered amplitude를 쓰면 t=0의 값이 0이므로 [Uₛ(+h)+Uₛ(−h)]/(2h²)다. U₂는 Taylor coefficient이며 두 번째 derivative 자체와 2배 차이가 있다는 점이 중요하다.

h=0.005와 0.01을 비교했다. 외부 scalar basis는 공통으로 유지하고 내부 index/basis와 scalar boundary coefficient만 바꿨다. 약한 대비의 scalar-only structure에 남아 있는 vector coefficient를 새로운 index의 Maxwell 해로 사용하지 않도록 guard를 두었다.

### 17.2 Detector diagnostic의 범위

계수 추출 후 detector 진단은 51.2 μm window, dx=0.1 μm, N=512, zdet=5 μm, NA=1.0에서 수행했다. 기존 forward algorithm과 맞추기 위해 forward-propagating spectrum으로 ψ₁,ψ₂를 만든 뒤 지수화하고 pupil을 적용했다.

따라서 이 ψ₂는 해당 propagating-spectrum field의 로그 계수다. Full Green-function xz 장에서 직접 계산한 Q와는 같은 계열의 진단이지만 동일한 데이터 경로는 아니다. Evanescent omission과 Rytov 절단을 완전히 따로 측정한 것은 아니다.

2차 계산은 experimental diagnostic이며, production born_rytov_voxel에 second-order Rytov 기능을 추가하지 않았다.

### 17.3 실제 대비 t=1에서의 결과

아래는 exact scalar reference에 대한 pupil 이후 complex scattered-field error이며 %다.

| θ | Born | 1차 Rytov | 2차 로그항 diagnostic | ‖ψ₂‖₂/‖ψ₁‖₂ |
|---:|---:|---:|---:|---:|
| 0° | 89.047 | 9.890 | 3.847 | 0.118669 |
| 5° | 89.273 | 9.957 | 3.878 | 0.119289 |
| 10° | 89.968 | 10.149 | 3.960 | 0.121368 |
| 15° | 91.139 | 10.479 | 4.138 | 0.125118 |

이 네 각도는 81-direction pupil 평균을 대신하지 않는다. 9.890%라는 θ=0 scalar-reference 값과 이전 θ=0 Maxwell co-reference의 9.869%도 reference가 다르다.

2차항을 추가한 뒤 약 10%가 약 4%가 되는 현상은 같은 detector 조건에서의 실질적인 모델 변경이다. 반면 10%에서 5%로 보였던 xz 그림의 비교는 앞 절의 metric 변경이다. 두 이야기를 혼동하지 않는다.

### 17.4 전체 대비 control

n(t)는 모든 θ에서 같다. t=0.1,0.25,0.5,1에 각각 n=1.338372877,1.342847397,1.350271980,1.365가 대응한다. t=0.5는 nₚ를 절반으로 만드는 것이 아니라 f를 절반으로 만드는 것이다.

| θ | t | 1차 Rytov 오차 (%) | 2차 diagnostic 오차 (%) |
|---:|---:|---:|---:|
| 0° | 0.10 | 0.8672 | 0.0339 |
| 0° | 0.25 | 2.1928 | 0.2140 |
| 0° | 0.50 | 4.5122 | 0.8798 |
| 0° | 1.00 | 9.8902 | 3.8472 |
| 5° | 0.10 | 0.8727 | 0.0342 |
| 5° | 0.25 | 2.2068 | 0.2158 |
| 5° | 0.50 | 4.5418 | 0.8872 |
| 5° | 1.00 | 9.9574 | 3.8776 |
| 10° | 0.10 | 0.8885 | 0.0349 |
| 10° | 0.25 | 2.2468 | 0.2208 |
| 10° | 0.50 | 4.6248 | 0.9069 |
| 10° | 1.00 | 10.1487 | 3.9600 |
| 15° | 0.10 | 0.9159 | 0.0364 |
| 15° | 0.25 | 2.3161 | 0.2301 |
| 15° | 0.50 | 4.7686 | 0.9458 |
| 15° | 1.00 | 10.4790 | 4.1382 |

[oblate_rytov_cause.csv](exp/oblate_rytov_cause.csv)에 Born과 ψ₂/ψ₁ norm ratio를 포함한 16개 전체 row를 보존했다.

### 17.5 원인에 대해 말할 수 있는 것

Scalar exact를 기준으로도 약 10.5%가 유지되고, independently extracted 2차 로그항이 오차를 줄이며, potential 대비를 낮추면 1차 오차가 줄어든다. 이는 로그장 expansion의 고차항 생략이 상당한 원인이라는 판단을 지지한다.

하지만 오차 norm이 10%에서 4%로 감소했다는 사실을 "전체 오차의 정확히 60%는 Q 때문"이라는 additive partition으로 바꾸지 않는다. 오차장들은 복소수로 간섭한다. 남은 약 4%를 특정 반사 차수, diffraction 또는 evanescent omission 하나에 전부 배정하지 않았다.

## 18. Q의 크기와 공간적 위치

### 18.1 단면 통계

[oblate_xz_nonlinear_term.csv](exp/oblate_xz_nonlinear_term.csv)는 exact scalar의 내부 xz 단면에서 계산한 ∣Q∣/∣f∣ 통계다.

| θ | Median | 90th percentile | 99th percentile | RMS | 내부 최소 ∣U∣ |
|---:|---:|---:|---:|---:|---:|
| 0° | 0.038445 | 0.196724 | 1.743499 | 0.478429 | 0.700020 |
| 5° | 0.038904 | 0.194322 | 1.685358 | 0.466157 | 0.708794 |
| 10° | 0.038548 | 0.206070 | 1.765847 | 0.495607 | 0.547076 |
| 15° | 0.037863 | 0.227345 | 1.714338 | 0.480821 | 0.612876 |

Q 통계에는 ∣U∣>0.1 조건을 준비했지만 모든 내부 점이 이를 만족해 실제 masking 비율은 0이었다. 이전 prolate의 exact field zero를 이번 oblate에서도 발견했다고 말하지 않는다.

이는 균일한 xz grid의 2차원 단면 통계이며 3차원 체적 평균이 아니다. RMS 약 0.5는 모든 지점에서 Q/f≈0.5라는 뜻도 아니다.

### 18.2 최대 위치: 꼭지점인가

[Q 위치 그림](exp/oblate_Q_xz_localization.png)은 내부만 표시하고 공통 log color scale 0.01–20을 사용한다. 외부에서는 f=0이므로 내부 비율의 연장처럼 색칠하지 않았다. Cyan +는 해당 저장 격자에서 가장 큰 값이다.

| θ | 최대 위치 x (μm) | 최대 위치 z (μm) | 최대 ∣Q∣/∣f∣ | ∣Q∣/∣f∣>1인 내부 점 비율 |
|---:|---:|---:|---:|---:|
| 0° | 0.00 | −0.20 | 14.3272 | 1.686% |
| 5° | 0.32 | −0.44 | 12.7342 | 1.727% |
| 10° | 0.72 | −0.28 | 14.4135 | 1.809% |
| 15° | 1.04 | −0.36 | 13.9679 | 1.904% |

이 위치들은 0.08 μm grid에서의 sampled maxima다. 연속 공간의 maximum을 최적화하거나 최고점 주변 grid를 더 세분한 결과가 아니다. 또한 y=0 단면 밖의 3차원 maximum을 조사한 것은 아니다.

가장 큰 영역은 z=±3 μm의 꼭지점이 아니라 내부 중앙 부근의 좁은 띠다. θ가 증가하면 그 띠가 +x 쪽으로 이동한다. 출사 쪽 +z 끝 근처에도 큰 secondary peak가 있지만, 해당 네 단면의 최고점은 내부 중앙에 있었다.

정상입사에서 ∣x∣<0.5 μm, ∣z∣<2 μm인 내부 영역은 단면의 Σ∣Q∣² 중 약 95.18%를 차지했다. 바깥 elliptic-radius shell 0.9<ρₑ<1은 약 1.73%였다. 이 비율은 Q의 제곱합 분포이지 detector error의 기여율이 아니다.

### 18.3 Q의 위치와 detector error의 위치는 같은 질문이 아니다

어느 지점에서 Q가 큰지는 exact field의 국소 성질이다. 그 지점이 detector에 만드는 complex error의 영향은 Green propagation, 다른 위치와의 위상 간섭, pupil에 따라 달라진다.

따라서 이 단계에서 내부 중앙부를 "detector 오차의 95%가 발생하는 곳"이라고 부르지 않는다. 그러한 인과 분할에는 Q의 spatial-source masking 또는 reflection-order field decomposition 같은 추가 계산이 필요하다.

## 19. 왜 내부 중앙에서 Q가 커지는가

### 19.1 먼저 exact derivative를 분해했다

Q 최대점의 log-amplitude gradient와 relative-phase gradient를 저장된 exact scalar 및 Cartesian gradient로 직접 계산했다. 정상입사 (x,z)=(0,−0.20) μm에서 다음 결과를 얻었다.

∣U∣ = 0.7041293921

∇α = (0,0,5.61290963) μm⁻¹

∇ΔΦ = (0,0,11.32853278) rad/μm

∣∇α∣²/f₀ = 2.82391705

−∣∇ΔΦ∣²/f₀ = −11.50331912

2i(∇α·∇ΔΦ)/f₀ = 11.39902084i

Q/f₀ = −8.67940207 + 11.39902084i

∣Q∣/f₀ = 14.32723618

따라서 실제 peak에서는 특히 z 방향 상대 위상 기울기와 amplitude–phase 교차항이 크다. 단순히 Q 그림의 밝은 부분이 전기장 intensity maximum이라는 뜻이 아니다. ∣U∣는 0에 가깝게 소실되지 않았지만 interference와 기울기 때문에 Q가 크다.

각 입사각의 Q peak에서 항별 분해는 다음과 같다. 모두 f₀로 나눈 값이며 마지막 열은 허수부 계수다.

| θ | ∣∇α∣²/f₀ | −∣∇ΔΦ∣²/f₀ | 2∇α·∇ΔΦ/f₀ |
|---:|---:|---:|---:|
| 0° | 2.823917 | −11.503319 | 11.399021 |
| 5° | 1.299710 | −11.478021 | 7.652614 |
| 10° | 0.362078 | −14.173441 | −4.122686 |
| 15° | 1.387222 | −12.594049 | −8.337228 |

### 19.2 Backward-wave 성분을 확인했다

정상입사 scalar field에서 ∣z∣<2.5 μm 구간에 Hann window를 적용하고 z 방향 FFT를 했다. 원래 dz=0.08 μm이며 FFT length 16,384로 zero padding했다. Zero padding은 spectrum을 더 촘촘하게 표시할 뿐 원래 공간 data의 독립 해상도를 늘리는 것은 아니다.

중앙축 x=0에서 양의 kz peak는 16.1212 μm⁻¹로 kₚ=16.1213307 μm⁻¹와 가깝고, 음의 kz peak는 −9.22785 μm⁻¹였다. 음/양 peak amplitude 비율은 약 0.11006이었다.

| x (μm), θ=0° | 음의 kz에 속하는 windowed 1D spectral squared norm 비율 |
|---:|---:|
| 0.00 | 0.064984 |
| 0.32 | 0.006553 |
| 0.64 | 0.002448 |
| 1.04 | 0.001103 |
| 2.00 | 0.000282 |
| 2.96 | 0.000132 |

같은 window를 pure forward plane wave에 적용한 음의 kz leakage 비율은 2.22e−8이었다. 따라서 관측된 음의 kz 성분은 단순 pure-plane-wave window leakage보다 훨씬 크고 중앙축 가까이 집중되어 있다는 증거다.

이것은 finite-window spatial-spectrum 진단이며, optical reflected power가 6.5%라는 측정이 아니다. 엄밀한 3차원 flux decomposition이나 내부 반사 차수 분해를 수행한 것도 아니다.

추가로 한 점의 U와 ∂zU를 ±kₚ 평면파로 algebraically 분해하는 간단한 probe도 계산했다. 이는 실제로 여러 각도의 파가 섞인 장에서 물리적 전진/반사 경로의 고유한 분해가 아니므로, 그 coefficient를 정확한 반사 amplitude나 Fresnel coefficient로 보고하지 않았다. 최종 해석의 근거는 전체 window의 spectrum과 아래 geometry 대조다.

### 19.3 실제 굴절률을 사용한 광선 추적

하부 경계에 +z로 입사하는 광선을 두고, Snell 법칙으로 내부 굴절 방향을 구한 다음 상부 경계에서 한 번 specular internal reflection시켰다. nₘ=1.335381534, nₚ=1.365, a=3 μm, b=5 μm를 그대로 사용했다.

입사 반경 1,2,3,4 μm의 대표 광선과 x→−x 대칭 광선을 추적했다. 굴절 방향의 unit norm, 다음 boundary 교점, 반사 후 축 교점이 다음 boundary 이전에 있는지 assertion으로 검사했다.

| 입사 반경 (μm) | 상부 반사점 x (μm) | 상부 반사점 z (μm) | 반사 후 축 교점 z (μm) |
|---:|---:|---:|---:|
| 1 | 0.984375 | 2.941286 | −1.038668 |
| 2 | 1.968731 | 2.757658 | −0.771842 |
| 3 | 2.953028 | 2.420881 | −0.240463 |
| 4 | 3.937078 | 1.849268 | 0.862068 |

상부 곡면은 내부에서 보면 오목한 반사면 역할을 하므로, 되돌아오는 광선이 내부 중앙축 근처로 모인다. 축대칭에서 같은 반경의 여러 방위각에서 온 광선이 축으로 모이는 현상은 내부 caustic 해석과 연결된다.

[반사 광선과 Q를 겹친 그림](exp/oblate_Q_reflected_ray_interpretation.png)에서 파란 선은 내부 투과 광선, 초록 선은 한 번 반사한 광선, cyan +는 exact Q의 sampled maximum이다. 광선이 모이는 z 구간이 Q가 큰 내부 중앙 영역과 겹친다.

처음의 빠른 geometry probe는 입구의 작은 굴절을 생략한 거의 평행한 내부 광선으로 축 교점을 추정했다. 최종 그림과 위 표는 그 단순화를 대체하여 Snell 굴절을 포함한 결과다.

### 19.4 해석의 강도와 남는 한계

현재 증거가 지지하는 해석은 **상부 곡면에서 생긴 내부 반사파의 축 부근 집중과 전진파와의 간섭이, 내부의 빠른 amplitude/phase variation과 큰 Q를 만든다**는 것이다.

Exact derivative에서 큰 기울기를 직접 확인했고, negative-kz 성분이 축 근처에 집중되며, 광선 geometry도 같은 내부 영역으로 반사광을 모은다. 단순히 그림 모양만 보고 꼭지점의 index jump 때문이라고 단정한 해석보다 구체적인 근거가 있다.

다만 ray calculation에는 Fresnel amplitude weight나 diffraction을 넣지 않았고, exact Maxwell/scalar field를 반사 차수별 Debye series로 분해하지 않았다. 따라서 초록 광선이 exact amplitude나 peak Q의 크기를 재현했다거나, 1회 반사가 10% detector error의 정해진 비율을 설명한다고 확정하지 않는다. Q maximum은 전진파와의 interference phase에 따라 geometric focus와 다른 위치에 있을 수 있다.

이 후속 분석은 [oblate_Q_mechanism.md](exp/oblate_Q_mechanism.md)에 별도로 보존했다. 내부반사 caustic에 관한 [관련 primary 연구](https://csuohio.elsevierpure.com/ws/portalfiles/portal/39955202/High-Order%20Interior%20Caustics%20Produced%20in%20Scattering%20of%20a%20Diagonal.pdf)는 일반적인 물리 배경이며, 다른 geometry의 논문으로 현재 parameter의 검증을 대신하지 않았다.

## 20. xz 장·미분·2차항 검증 결과

| 검사 | 상대 차이 또는 normalized residual |
|---|---:|
| Near vector, (109,55) 대 (114,59) | 2.468353e−12 |
| Near scalar, 같은 cutoff 비교 | 2.721561e−12 |
| Potential derivative Born 대 기존 surface-integral Born | 1.682279e−5 |
| Born의 h=0.005 대 0.01 | 4.732921e−5 |
| 2차 coefficient의 h refinement | 2.822297e−5 |
| Batched Ex 대 public spheroid_eval | 7.343261e−16 |
| Scalar Cartesian gradient 대 finite difference | 1.613801e−8 |
| 정확한 logarithmic Helmholtz PDE residual | 5.428536e−7 |

Near cutoff는 257개 grid point에서 비교했다. Public evaluator 및 derivative/PDE 검사는 내부·외부 8개 point와 네 입사각에서 수행했다. Derivative finite difference step은 2e−5 μm다. 모든 35,376개 점에서 개별 PDE residual을 검사한 것은 아니다.

Independent Born surface-integral 비교는 별도의 내부·외부 8개 점에서 기존 함수를 tol=1e−7로 실행하여 수행했다. 약한 potential derivative의 finite h 차이는 10⁻⁵ 수준으로 확인되었으며, 2.8 같은 field magnitude나 약 10% 근사 오차와 비교해 작다.

완료 증거는 다음과 같다.

- oblate_xz.log: test_sph_coords: PASS, test_coords_vecwave: PASS, COORDINATE_CHECKS_PASS.
- 같은 log: XZ_MAPS_SAVED, 최종 OBLATE_XZ_AND_CAUSE_PASS.
- oblate_xz_derivative_validation.log: XZ_PUBLIC_AND_PDE_PASS.
- oblate_final_checks.log: OBLATE_IMPLEMENTATION_CHECKS_PASS, OBLATE_FINAL_CHECKS_PASS.
- 각 검증 수치: [xz validation CSV](exp/oblate_xz_validation.csv), [derivative validation CSV](exp/oblate_xz_derivative_validation.csv).

추가 Python artifact check는 exact/Born/Rytov의 12개 field panel이 finite이고 176×201×4에 맞는지, θ가 0,5,10,15인지, contrast CSV 16행과 plot-error CSV 24행이 있는지, PNG가 정상적으로 열리고 PDF header가 유효한지 확인했다. 16개 contrast 조건 모두에서 이번 2차 diagnostic의 오차가 1차보다 작았다.

코드 검토에서는 exterior common-basis 차분과 interior direct-field 차분, cosθ projection, Cartesian gradient metric factor, U₂의 2h² 분모, MATLAB column-major reshape, 공유 색상척도, scalar-only guard를 확인했다. 이 read-only 검토와 실제 numerical execution의 PASS는 다른 종류의 증거이며 함께 기록한다.

## 21. 실행 환경

MATLAB은 /Applications/MATLAB_R2024a.app/bin/matlab의 R2024a Update 9였다.

계산은 -nojvm -nodesktop -nosplash -r로 시작하고, session-local restoredefaultpath 후 필요한 project path를 넣었다. 최종 판단은 try/catch, exit code, assertion, 완료 marker를 사용했다.

전역 MATLAB 설정이나 다른 사용자 process는 변경하지 않았다. 이전 prolate 및 scalar-pilot 기록에 등장하는 중단된 public solve는 완료된 batch 검증과 구분해 보존한다.

수치 계산을 MATLAB에서 완료한 뒤 저장된 MAT를 Python으로 그렸다. 최종 PNG/PDF 생성은 완료되었다.

Plot 전용 임시 가상환경을 /tmp/mie-oblate-plot-env에 만들고 NumPy 2.4.6, SciPy 1.17.1, Matplotlib 3.11.2를 사용했다. 시스템 Python 환경의 library 설치를 대신한 task-local 환경이며, /tmp 환경이 장기간 보존된다는 보장은 없다. 코드나 PNG/MAT를 재사용하려면 해당 library가 있는 Python 환경을 준비하면 된다.

가장 무거운 작업은 exact 및 약한 대비의 interior near-field sum이었다. 공통 reference batch는 수십 초에 끝났지만, 35,376점의 공간장과 약한 대비 ±h 내부장을 평가하는 단계가 더 오래 걸렸다. 실행시간은 현재 컴퓨터와 MATLAB 상태에 의존하며 algorithm complexity의 보편적 benchmark로 사용하지 않는다.

중간의 plotting probe나 Python bytecode처럼 이번 작업에서 만든 불필요한 임시 파일은 정리했다. 사용자 데이터나 다른 project tree를 정리 대상으로 삼지 않았다.

## 22. 보존 파일과 각 파일의 역할

### 22.1 Production·기존 구현 및 회귀검사

| 파일 | 역할 |
|---|---|
| [plan.md](plan.md) | 연속 spheroid Born/Rytov의 수학·API·수렴 명세 |
| [born_rytov_spheroid.m](born_rytov_spheroid.m) | 임의 내부·외부 점의 full-Green scalar surface-integral forward |
| [born_rytov_voxel.m](born_rytov_voxel.m) | FFT/Ewald 및 올바른 pupil 순서의 detector forward |
| [compare_with_spheroid.m](compare_with_spheroid.m) | 기존 reference 비교 도구 |
| [test_born_rytov.m](tests/test_born_rytov.m) | Sphere series, volume quadrature, far field, Rytov 및 failure-policy 검사 |
| [test_born_rytov_voxel.m](tests/test_born_rytov_voxel.m) | Voxel/Ewald 및 NA 순서 회귀검사 |
| [run_all.m](tests/run_all.m) | 현재는 test_born_rytov를 호출하는 runner |
| [analyze_rytov_spheroid.m](tests/analyze_rytov_spheroid.m) | 이전 prolate detector comparison |

현재 forward/tests/run_all.m이 voxel 및 모든 experimental stage까지 포함하는 통합 runner라고 생각하면 안 된다. 필요한 함수를 명시적으로 따로 실행한다.

### 22.2 이전 prolate experiment

| 파일 묶음 | 내용 |
|---|---|
| [scalar_prolate.m](exp/scalar_prolate.m), [scalar_prolate_eval.m](exp/scalar_prolate_eval.m), [scalar_far_amplitude.m](exp/scalar_far_amplitude.m) | 기존 axial prolate scalar reference |
| [test_scalar_prolate.m](exp/test_scalar_prolate.m) | Zero contrast, sphere 등 scalar reference 검사 |
| [diagnose_rytov.m](exp/diagnose_rytov.m) | Contrast/distance spectrum 및 Rytov order 분석 |
| [diagnose_rytov_term.m](exp/diagnose_rytov_term.m) | Exact logarithmic gradient 및 Q 분석 |
| [diagnose_rytov_geometry.m](exp/diagnose_rytov_geometry.m) | Prolate/sphere/slab control |
| [locate_field_zero.m](exp/locate_field_zero.m) | Scalar total-field zero 탐색 |
| [report.md](exp/report.md) | 이전 prolate 원인 보고서 |
| [verification.log](exp/verification.log), [diagnosis.log](exp/diagnosis.log) | 완료 검증과 개발 중 실패의 보존 기록 |
| [contrast_distance_N512.csv](exp/contrast_distance_N512.csv), [contrast_distance_N1024.csv](exp/contrast_distance_N1024.csv) | 거리·대비·Rytov 차수별 결과 |
| [geometry_control.csv](exp/geometry_control.csv), [field_zero.csv](exp/field_zero.csv) | Geometry 대조 및 zero의 좌표 |
| [omitted_term_dx0.1.csv](exp/omitted_term_dx0.1.csv), [omitted_term_dx0.05.csv](exp/omitted_term_dx0.05.csv), [omitted_term_dx0.025.csv](exp/omitted_term_dx0.025.csv) | 이전 Q grid refinement |

이전 spectral_N512/N1024.mat, term_t01/t10_dx*.mat, slab_control.mat도 exp에 남아 있다. 이들 prolate MAT와 이름이 oblate_로 시작하는 새 MAT는 같은 parameter의 checkpoint가 아니다.

### 22.3 Oblate sweep와 검증

| 파일 | 내용 |
|---|---|
| [oblate_na_plan.md](exp/oblate_na_plan.md) | 승인된 detector-sweep 범위와 완료 목록 |
| [oblate_na_sweep.m](exp/oblate_na_sweep.m) | check/reference/sweep/fieldcheck/report/xz/xz_validate의 재현 script |
| [oblate_report.md](exp/oblate_report.md) | 81-direction sweep의 결과·metric·수치 한계 |
| [oblate_angles.csv](exp/oblate_angles.csv) | 264개의 sweep/control row |
| [oblate_summary.csv](exp/oblate_summary.csv) | 33개의 configuration/metric summary |
| [oblate_baseline.csv](exp/oblate_baseline.csv) | 81-direction baseline |
| [oblate_window_51.csv](exp/oblate_window_51.csv) | 81-direction 51.2 μm window |
| [oblate_sampling_0p05.csv](exp/oblate_sampling_0p05.csv) | 81-direction fine sampling |
| [oblate_padding_3.csv](exp/oblate_padding_3.csv) | 7-direction padding control |
| [oblate_voxel_0p0667.csv](exp/oblate_voxel_0p0667.csv) | 7-direction voxel pitch control |
| [oblate_voxel_window_51.csv](exp/oblate_voxel_window_51.csv) | 7-direction voxel FOV control |
| [oblate_reference_validation.csv](exp/oblate_reference_validation.csv) | Boundary/cutoff/far-limit 검증 |
| [oblate_field_convergence.csv](exp/oblate_field_convergence.csv) | 동일 물리 주파수에서 field refinement |
| [oblate_angular_convergence.csv](exp/oblate_angular_convergence.csv) | Angular 평균 coarsening 검사 |
| [oblate_reference_0.mat](exp/oblate_reference_0.mat), [oblate_reference_1.mat](exp/oblate_reference_1.mat) | Coarse/fine compact far-field reference |
| [oblate_checks.log](exp/oblate_checks.log), [oblate_reference.log](exp/oblate_reference.log) | 구현 및 reference 실행 |
| [oblate_sweep.log](exp/oblate_sweep.log), [oblate_field_convergence.log](exp/oblate_field_convergence.log) | Sweep와 수치 refinement 실행 |
| [oblate_final_checks.log](exp/oblate_final_checks.log) | 최종 집계 및 assertion 완료 |

### 22.4 xz 그림과 Q 해석

| 파일 | 내용 |
|---|---|
| [oblate_xz_report.md](exp/oblate_xz_report.md) | 전체 Ex 그림, Rytov 원인 및 검증 보고서 |
| [plot_oblate_xz.py](exp/plot_oblate_xz.py) | 4×3 그림과 plot-error CSV 생성 |
| [oblate_Ex_xz_4x3.png](exp/oblate_Ex_xz_4x3.png), [PDF](exp/oblate_Ex_xz_4x3.pdf) | 요청된 12개 field panel |
| [oblate_xz_fields.mat](exp/oblate_xz_fields.mat) | Exact/Born/Rytov 복소 Ex, Q, scalar/incident/Born field 및 grid |
| [oblate_xz_exact.mat](exp/oblate_xz_exact.mat) | Exact Ex, scalar U, Cartesian gradient, Q 및 grid |
| [oblate_xz_plot_errors.csv](exp/oblate_xz_plot_errors.csv) | 전체/내부 영역의 complex 및 magnitude error 24행 |
| [oblate_xz_metrics.csv](exp/oblate_xz_metrics.csv) | 네 각도의 Q RMS, minimum U, total Ex error |
| [oblate_xz_nonlinear_term.csv](exp/oblate_xz_nonlinear_term.csv) | Q median/percentile/RMS |
| [oblate_rytov_cause.csv](exp/oblate_rytov_cause.csv) | 16 angle/contrast combination |
| [oblate_contrast_t0.1.mat](exp/oblate_contrast_t0.1.mat), [t0.25](exp/oblate_contrast_t0.25.mat), [t0.5](exp/oblate_contrast_t0.5.mat) | 각 대비의 네 입사각 exact scalar detector spectrum |
| [oblate_xz_validation.csv](exp/oblate_xz_validation.csv), [derivative validation](exp/oblate_xz_derivative_validation.csv) | 수치 검사 |
| [oblate_xz.log](exp/oblate_xz.log), [derivative log](exp/oblate_xz_derivative_validation.log) | 성공한 xz/미분 실행 |
| [oblate_Q_xz_localization.png](exp/oblate_Q_xz_localization.png) | 네 입사각의 내부 Q 위치 colormap |
| [oblate_Q_reflected_ray_interpretation.png](exp/oblate_Q_reflected_ray_interpretation.png) | 정상입사 Q와 굴절·반사 광선 대조 |
| [oblate_Q_mechanism.md](exp/oblate_Q_mechanism.md) | Gradient/FFT/ray 기반 물리적 해석 |
| [history.md](history.md) | 이 대화의 통합 이력 |

Q localization과 ray-overlay 그림 및 nonlinear-term percentile CSV는 저장된 MAT를 읽는 단발 Python 분석으로 만들었다. 독립 plot_Q.py를 추가하지 않았으며, plot_oblate_xz.py가 이 두 Q 그림까지 생성한다고 설명하면 안 된다. 원시 MAT, 위의 plotting 설정, 아래 수치 재현 코드와 그림 자체가 남아 있다.

## 23. 재현 절차

### 23.1 핵심 MATLAB stage

저장소 root에서 MATLAB에 다음 순서로 명령을 준다.

```matlab
restoredefaultpath;
addpath('forward/exp');
oblate_na_sweep('check');
oblate_na_sweep('reference');
oblate_na_sweep('sweep');
oblate_na_sweep('fieldcheck');
oblate_na_sweep('report');
oblate_na_sweep('xz');
oblate_na_sweep('xz_validate');
```

Reference stage는 parameter가 일치하는 compact checkpoint가 있으면 재사용한다. Sweep/report stage는 그 checkpoint와 CSV를 전제로 한다. xz 및 xz_validate는 물리 parameter로 full batch reference를 다시 만든다. 코드를 변경한 뒤 cached 결과가 새 구현을 검증한다고 자동으로 생각하면 안 된다.

현재 host에서 사용한 실행 형태는 다음과 같다. 이 명령은 shell 기준으로 실행하며 해당 log 파일에 출력한다.

```sh
/Applications/MATLAB_R2024a.app/bin/matlab -nojvm -nodesktop -nosplash -r "restoredefaultpath; addpath('forward/exp'); try, oblate_na_sweep('xz'); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0)" > forward/exp/oblate_xz.log 2>&1
```

좌표 fix의 focused regression은 spheroid-analytic-forward와 그 tests를 MATLAB path에 넣은 뒤 test_sph_coords, test_coords_vecwave를 실행한다. Voxel operator 검사는 forward/tests의 test_born_rytov_voxel을 별도로 실행한다.

### 23.2 그림 재생성

NumPy, SciPy, Matplotlib가 있는 Python 환경에서 다음을 실행한다.

```sh
python forward/exp/plot_oblate_xz.py
```

이번 실행에서는 python 대신 /tmp/mie-oblate-plot-env/bin/python을 사용했다. Script는 exp/oblate_xz_fields.mat을 읽어 PNG, PDF와 24-row plot-error CSV를 생성한다.

MATLAB 배열의 spatial flattening은 column-major다. Grid 복원은 항상 다음 규칙을 따른다.

field[:, angle].reshape((len(z), len(x)), order="F")

Scalar Ex는 scalar × cosθ이고, Ex 배열에 이미 들어간 cosθ를 다시 곱하지 않는다.

### 23.3 Q peak 및 항 분해 재현

다음 코드는 새 MATLAB solve 없이 저장된 exact 장으로 최대 위치와 항별 값을 재현한다. 모든 단위는 저장된 μm grid 기준이다.

```python
from pathlib import Path
import numpy as np
from scipy.io import loadmat

d = loadmat(Path("forward/exp/oblate_xz_exact.mat"), squeeze_me=True)
X, Z = np.meshgrid(d["x"], d["z"])
x, z = X.ravel(order="F"), Z.ravel(order="F")
inside = d["inside"].astype(bool)
km = 2 * np.pi * 1.335381534 / 0.532

for j, angle in enumerate(d["theta_deg"]):
    theta = np.deg2rad(angle)
    ki = km * np.array([np.sin(theta), 0, np.cos(theta)])
    u = d["scalar"][:, j]
    g = d["grad"][:, :, j] / u[:, None] - 1j * ki
    q = np.sum(g * g, axis=1)
    assert np.allclose(q, d["Q"][:, j])
    index = np.argmax(np.where(inside, np.abs(q), -1))
    ga, gp = g[index].real, g[index].imag
    print(angle, x[index], z[index], abs(q[index]) / d["f0"])
    print("gradient terms / f:",
          np.dot(ga, ga) / d["f0"],
          -np.dot(gp, gp) / d["f0"],
          2 * np.dot(ga, gp) / d["f0"])
```

### 23.4 Backward-spectrum 진단 재현

위와 같은 d를 읽은 상태에서 실행한다. 결과는 finite-window spatial spectrum이며 power reflection coefficient가 아니다.

```python
z, x = d["z"], d["x"]
dz = z[1] - z[0]
kp = 2 * np.pi * 1.365 / 0.532
field = d["scalar"][:, 0].reshape((len(z), len(x)), order="F")
window_points = abs(z) < 2.5
window = np.hanning(window_points.sum())
kz = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(16384, dz))

for target_x in [0, 0.32, 0.64, 1.04, 2, 3]:
    ix = np.argmin(abs(x - target_x))
    signal = field[window_points, ix] * window
    spectrum = np.fft.fftshift(np.fft.fft(signal, 16384))
    fraction = np.sum(abs(spectrum[kz < 0])**2) / np.sum(abs(spectrum)**2)
    print(x[ix], fraction)

pure = np.exp(1j * kp * z[window_points]) * window
pure_spectrum = np.fft.fftshift(np.fft.fft(pure, 16384))
print("pure-wave negative-kz leakage:",
      np.sum(abs(pure_spectrum[kz < 0])**2) / np.sum(abs(pure_spectrum)**2))
```

### 23.5 굴절·반사 광선의 교점 재현

이것은 phase/amplitude를 계산하는 wave solver가 아니라 geometry control이다.

```python
import numpy as np

a, b, nm, nparticle = 3.0, 5.0, 1.335381534, 1.365
metric = np.array([1 / b**2, 1 / a**2])

def next_boundary(point, direction):
    distance = -2 * np.dot(point * metric, direction) / np.dot(direction * metric, direction)
    assert distance > 0
    return point + distance * direction, distance

for radius in [1, 2, 3, 4]:
    entry = np.array([radius, -a * np.sqrt(1 - radius**2 / b**2)])
    inward = -entry * metric
    inward /= np.linalg.norm(inward)
    incident = np.array([0.0, 1.0])
    tangent = (nm / nparticle) * (incident - np.dot(incident, inward) * inward)
    transmitted = tangent + np.sqrt(1 - np.dot(tangent, tangent)) * inward
    assert abs(np.linalg.norm(transmitted) - 1) < 1e-12
    upper, _ = next_boundary(entry, transmitted)
    outward = upper * metric
    outward /= np.linalg.norm(outward)
    reflected = transmitted - 2 * np.dot(transmitted, outward) * outward
    end, length = next_boundary(upper, reflected)
    distance = -upper[0] / reflected[0]
    assert 0 < distance < length
    crossing = upper + distance * reflected
    print(radius, upper[0], upper[1], crossing[1])
```

## 24. 이 작업에서 바뀐 파일과 repository 상태

이번 oblate 작업에서 핵심 실험 코드는 exp/oblate_na_sweep.m, plotting 코드는 exp/plot_oblate_xz.py에 들어갔다. CSV/MAT/PNG/PDF/report/log를 함께 보존했다. 이어 사용자의 명시적 요청으로 forward/history.md를 새로 만들었다.

Raw near-field 계산을 위해 shared solver의 sph_coords.m 한 줄과 test_sph_coords.m의 regression을 수정했다. n, geometry, 편광 convention을 바꿔 오류를 피한 것은 아니다. Production born_rytov_voxel에 2차항을 추가하거나 inverse 알고리즘을 새 모델에 맞춰 수정하지 않았다.

이 문서 작성 시 확인한 branch는 main, HEAD는 51094f1 ("docs(spheroid): record scalar validation scope")였다. forward/와 inverse/ 등은 untracked였고 spheroid-analytic-forward/history.md에는 이미 별도 변경이 있었다. 이전 sibling history 안의 과거 HEAD/commit 상태는 그 기록 시점의 정보이며 현재 상태와 구분한다.

이번 history 작성 요청에서는 sibling history를 덮거나 append하지 않고 forward/history.md만 생성·편집했다. 이 문서 작성 과정에서 commit, push, merge, global configuration 변경, memory database 갱신은 하지 않았다.

## 25. 확정된 결론과 아직 하지 않은 일

### 25.1 현재 증거로 확정할 수 있는 것

1. 지정한 oblate 형상과 굴절률, 기존 pol=[1,0] convention, illumination NA≤0.5와 detector NA=1.0에서 Born/Rytov 성능을 재현 가능한 81-direction 계산으로 평가했다.
2. 연속형상 Rytov의 Maxwell co-relative complex scattered-field weighted mean은 약 10.509%이고, scalar exact와 비교해도 약 10.525%였다.
3. 같은 detector에서 scalar exact의 co 오차는 약 0.705%, full-vector 오차는 약 3.117%였다. Scalar reduction과 Rytov 절단은 다른 근사다.
4. Baseline voxel의 더 작은 9.751%에는 수치 오차의 상쇄가 포함되어 있으며, 그 숫자만으로 물리 근사가 더 좋다고 해석하지 않는다.
5. 요청된 네 각도의 raw total Ex를 4×3 그림으로 저장했다.
6. 같은 xz 배열에서 복소장 오차는 약 8–9%, 크기만의 오차는 약 5–6%였으며, 10%→5%를 알고리즘 개선으로 해석한 혼동을 정정했다.
7. Exact scalar의 Q는 대부분 작은 값이지만 내부 중앙의 좁은 영역에서 커졌고, 해당 sampled maxima는 꼭지점이 아니었다.
8. 2차 로그항과 대비 control은 고차항 절단이 상당한 오차원이라는 판단을 지지했다.
9. Negative-kz 및 광선 추적 결과는 내부 반사파의 집중과 interference라는 물리적 해석을 지지했다.
10. 구현·수치 검증과 그림·원자료·완료 log를 보존했다.

### 25.2 아직 완료하거나 증명하지 않은 것

- 실제 biological cell, 세포질/핵/membrane의 heterogeneous 또는 layered 모델 검증.
- Scalar/Born/Rytov의 모든 wavelength·size·contrast·aspect ratio·orientation·absorption 범위.
- 연속된 모든 illumination angle의 엄밀한 최대 오차 bound.
- 실제 objective의 polarization transport, aberration, analyzer 및 detector response.
- Evanescent omission과 Rytov order truncation의 완전한 독립 오차 분할.
- Q의 공간별 source를 detector까지 전파해 각 위치의 error contribution을 계산하는 분석.
- 반사 차수별 exact field decomposition, Fresnel-weighted wave caustic 재현.
- 3차원 전체 Q map의 maximum 또는 현재 0.08 μm grid peak의 별도 spatial refinement.
- 2차 Rytov의 production 구현 및 모든 조건에서의 안정성 보장.
- 현재 voxel field가 1% 이내로 수렴했다는 증명.
- 새 oblate 모델을 사용한 inverse reconstruction 및 reconstruction-error 평가.
- 모든 기존 test suite를 현재 tree에서 한 번에 다시 돌린 통합 결과.

기존 sibling 기록에는 Rytov pupil 수정 후 inverse test/demo의 sampling 및 logarithm convention이 맞지 않는 문제가 남아 있다고 되어 있다. 이번 작업은 inverse를 수정하거나 재검증하지 않았으므로 그 항목을 해결 완료로 표시하지 않는다.

이전의 72-case comparison sweep나 537-case 전체 solver stress campaign도 이 oblate 작업에서 새로 완료하지 않았다. 일부 focused test와 experimental assertion의 성공을 전체 campaign 완료로 확대하지 않는다.

## 26. 사용한 이론 자료와 출처의 한계

- [Müller, Schürmann, Guck, The Theory of Diffraction Tomography, arXiv:1507.00466v3](https://arxiv.org/pdf/1507.00466): scalar Helmholtz, Born/Rytov, 특히 §3.2의 로그장 방정식과 생략항. [로컬 PDF](../ref/background/1507.00466v3.pdf).
- [Kak and Slaney, Chapter 6: Tomographic Imaging with Diffracting Sources](../ref/background/SIAMB0000017_chapter-10_Chapter_6_Tomographic_Imaging_with_Diffracting_Sources.pdf): Green 함수, volume integral, Born/Rytov 및 diffraction tomography의 이론 배경.
- [spheroid-analytic-forward/plan.md](../spheroid-analytic-forward/plan.md), [그 개발 이력](../spheroid-analytic-forward/history.md): geometry·편광·시간 convention, Maxwell reference 및 특수함수 검증의 배경.
- [Tomocube/Rytov reference 1](../ref/tomocube-rytov/oe-23-13-16933.pdf), [reference 2](../ref/tomocube-rytov/s41598-018-25886-8.pdf): 사용자가 제시한 biological/ODT 문헌 배경. 현재 축 길이와 모든 index를 그대로 뽑아 온 특정 표를 이 대화에서 새로 인증한 것은 아니다.
- [High-order interior caustics produced in scattering of a diagonally incident plane wave by a circular cylinder](https://csuohio.elsevierpure.com/ws/portalfiles/portal/39955202/High-Order%20Interior%20Caustics%20Produced%20in%20Scattering%20of%20a%20Diagonal.pdf): 내부반사와 caustic의 일반 배경. 대상 geometry가 다르므로 현재 oblate의 정량 검증 자료로 사용하지 않는다.

기존 forward/plan.md에는 reference 위치를 forward/ref로 부르는 과거 경로 설명이 있으나, 이 문서에서 연결한 실제 background PDF 위치는 repository의 ref/background다. 원본 plan을 이 history 요청 중에 수정하지는 않았다.

이 파일의 수치 결론은 연결된 코드·CSV·MAT·log에 근거한다. 문헌의 일반 이론, 과거 console-only 기록, 현재 재현 가능한 실험, 후속 물리적 추론을 같은 수준의 증거로 취급하지 않는다.
