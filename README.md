# fmriroi

Layered fMRI ROI Clustering Framework

---

## 개요 (Overview)

`fmriroi`는 ROI 기반 fMRI 분석을 위한 모듈형 Python 패키지입니다.

본 패키지는 다음과 같은 **5단계 레이어 아키텍처**를 기반으로 설계되었습니다:

1. **Data Model**
2. **Transform (전처리 / 차원축소)**
3. **Similarity Computation (유사도 계산)**
4. **Clustering**
5. **Visualization & Reporting**

데이터는 아래와 같은 단방향 흐름을 따릅니다:

```

FmriVoxelData
→ PreprocessedData
→ SimilarityResult
→ ClusteringResult

````

각 레이어는 하나의 클래스 객체로 표현되며,  
실제 계산 알고리즘은 별도의 모듈로 분리되어 있습니다.

---

## 설계 철학 (Design Philosophy)

- 상태 객체와 계산 알고리즘의 명확한 분리
- 단방향 데이터 흐름 구조
- 각 단계별 메타데이터 자동 추적
- 연구 확장성을 고려한 모듈형 설계
- Jupyter Notebook 및 Gradio GUI 확장 가능

---

## 설치 방법 (개발 모드)

```bash
git clone https://github.com/<your-username>/fmriroi.git
cd fmriroi
pip install -e .
````

---

## 사용 예시 (Quick Example)

```python
import numpy as np
from fmriroi import FmriVoxelData

# 예시 데이터
X = np.random.randn(200, 500)      # (T x V)
coords = np.random.randint(0, 64, size=(500, 3))

fmri = FmriVoxelData(X, coords).validate()

prep = fmri.preprocess("bspline", n_basis=10)

sim = prep.compute_similarity("correlation", to_affinity=True)

clu = sim.cluster("spectral_craddock26", n_clusters=100)

clu.plot_voxels_3d()

score = clu.report_silhouette()
```

---

# 아키텍처 설명

## Layer 1 – Data Model

입력 데이터를 내부 표준 포맷으로 정리하고 검증합니다.

주요 기능:

* `make_fmri_tensor`
* `make_roi_table`
* `validate_data`
* `split_train_test`

주요 클래스:

* `FmriVoxelData`

---

## Layer 2 – Transform (전처리)

스무딩 및 차원축소를 수행합니다.

지원 기능:

* Gaussian smoothing
* B-spline basis expansion
* Fourier basis expansion (예정)
* fPCA
* fICA (예정)
* Standardization

주요 클래스:

* `PreprocessedData`

---

## Layer 3 – Similarity Computation

특징 벡터 또는 시계열 간 유사도/거리 행렬을 생성합니다.

지원 방법:

* Pearson correlation
* Partial correlation
* Cosine similarity
* L2 distance
* SRVF distance (예정)
* Affinity 변환 (RBF 등)

주요 클래스:

* `SimilarityResult`

---

## Layer 4 – Clustering

유사도 행렬 또는 특징 공간을 기반으로 클러스터링을 수행합니다.

지원 방법:

* Spectral clustering (Craddock 26-neighbor 제약)
* Spectral clustering (비제약)
* K-means (예정)
* Hierarchical clustering
* Cluster 후처리

주요 클래스:

* `ClusteringResult`

---

## Layer 5 – Visualization & Reporting

시각화 기능:

* Similarity heatmap
* Dendrogram
* 3D voxel scatter plot
* Cluster 색상 시각화
* Figure 저장 기능

평가 지표:

* Silhouette score
* DICE coefficient

---

# 프로젝트 구조

```

src/fmriroi/
├── data/
├── transform/
├── similarity/
├── cluster/
├── viz/
├── metrics/
├── utils/

```

* `data/`는 레이어 상태 객체를 포함합니다.
* 다른 모듈은 실제 계산 엔진을 포함합니다.
* `__init__.py`는 외부에 노출되는 API를 정리합니다.

---

# 개발 가이드

권장 개발 흐름:

1. 계산 알고리즘은 각 레이어 모듈에 구현
2. 레이어 클래스 메서드를 통해 연결
3. 기능 단위로 커밋
4. 새로운 알고리즘은 feature 브랜치에서 개발

---

# 향후 계획 (Roadmap)

* Fourier basis 확장
* fICA 구현
* SRVF 최적화
* GPU 가속 지원
* Gradio 기반 GUI 통합

---

# 라이선스

(To be specified)

---

# 작성자

Hyunseok Yoon

```
