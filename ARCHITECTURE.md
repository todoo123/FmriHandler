# Architecture

## TL;DR (핵심 요약)

- 데이터 흐름:
  `FmriVoxelData → PreprocessedData → SimilarityResult → ClusteringResult`
- 각 레이어 클래스는 "흐름 제어(오케스트레이션)"만 담당한다.
- 실제 계산은 별도 모듈(엔진 함수)에서 수행한다.
- 순환 import를 허용하지 않는다.
- 데이터 흐름은 단방향으로만 진행된다.

---

# 1. 전체 구조 개요

본 프로젝트는 다음과 같은 5단계 레이어 아키텍처를 따른다.

1. Data Model
2. Transform (전처리 / 차원 축소)
3. Similarity Computation (유사도 계산)
4. Clustering
5. Visualization & Reporting

데이터는 아래와 같은 단방향 파이프라인을 따른다.

FmriVoxelData
→ PreprocessedData
→ SimilarityResult
→ ClusteringResult


각 단계는 새로운 객체를 반환하며, 이전 단계로 되돌아가지 않는다.

---

# 2. 레이어별 책임 분리

## 2.1 data/

상태 객체(레이어 클래스)를 정의하는 모듈.

클래스:
- `FmriVoxelData`
- `PreprocessedData`
- `SimilarityResult`
- `ClusteringResult`

역할:
- 데이터 및 메타데이터 보관
- 다음 단계로의 디스패치
- 엔진 함수 호출
- 다음 레이어 객체 생성

금지 사항:
- 무거운 수치 계산 로직 구현
- 알고리즘 세부 구현 포함

---

## 2.2 transform/

전처리 및 차원 축소 엔진 모듈.

예시:
- Gaussian smoothing
- B-spline basis expansion
- fPCA
- Standardization

규칙:
- data 모듈을 import하지 않는다.
- numpy 배열 기반 순수 함수로 구현한다.
- 메타데이터를 관리하지 않는다.

---

## 2.3 similarity/

유사도 및 거리 계산 엔진 모듈.

예시:
- Pearson correlation
- Partial correlation
- Cosine similarity
- L2 distance
- Affinity 변환

규칙:
- 순수 함수만 작성
- clustering 모듈에 의존하지 않는다.

---

## 2.4 cluster/

클러스터링 알고리즘 모듈.

예시:
- Spectral clustering (Craddock 26-neighbor 제약)
- Spectral clustering (비제약)
- Hierarchical clustering
- K-means (예정)
- Postprocessing

규칙:
- 입력은 similarity/affinity 행렬 또는 특징 행렬
- 출력은 cluster label 배열

---

## 2.5 viz/ 및 metrics/

시각화 및 평가 지표 모듈.

예시:
- Heatmap
- Dendrogram
- 3D voxel plot
- Silhouette score
- DICE metric

규칙:
- 상태 객체를 수정하지 않는다.
- 계산 모듈과 순환 의존을 만들지 않는다.

---

# 3. 의존성 규칙 (중요)

허용:

- `data/` → transform/similarity/cluster/viz/metrics import 가능
- transform/similarity/cluster → utils import 가능

금지:

- transform/similarity/cluster가 data를 import
- 레이어 간 상향 의존
- 순환 import

---

# 4. Method Dispatch 패턴

각 레이어 클래스는 얇은 오케스트레이션 메서드를 가진다.

- `FmriVoxelData.preprocess(method, **kwargs)`
- `PreprocessedData.compute_similarity(method, **kwargs)`
- `SimilarityResult.cluster(method, **kwargs)`
- `ClusteringResult.plot_*() / report_*()`

디스패치는 내부 매핑 딕셔너리로 관리한다.

예시:

```python
_PREPROCESS_MAP = {
    "gaussian": gaussian_kernel_smooth,
    "bspline": basis_expand_bspline,
}

