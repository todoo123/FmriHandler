# fmriroi

Layered fMRI ROI Clustering Framework

---

## Overview

`fmriroi` is a modular and extensible Python package for ROI-based fMRI analysis.

It implements a **5-layer pipeline architecture**:

1. **Data Model**
2. **Transform (Preprocessing / Dimension Reduction)**
3. **Similarity Computation**
4. **Clustering**
5. **Visualization & Reporting**

The design enforces a clean, unidirectional workflow:

```

FMRIData
→ PreprocessedData
→ SimilarityResult
→ ClusteringResult

````

Each layer is represented by a dedicated class, while computational engines are implemented in separate modules for clarity and extensibility.

---

## Design Philosophy

- Clear separation between **state objects** and **algorithm engines**
- Layered, one-directional data flow
- Minimal but explicit metadata tracking
- Research-oriented modular design
- Suitable for Jupyter + GUI (Gradio) integration

---

## Installation (Development Mode)

Clone the repository and install in editable mode:

```bash
git clone https://github.com/<your-username>/fmriroi.git
cd fmriroi
pip install -e .
````

---

## Quick Example

```python
import numpy as np
from fmriroi import FMRIData

# Example data
X = np.random.randn(200, 500)        # (T x V)
coords = np.random.randint(0, 64, size=(500, 3))

fmri = FMRIData(X, coords).validate()

prep = fmri.preprocess("bspline", n_basis=10)

sim = prep.compute_similarity("correlation", to_affinity=True)

clu = sim.cluster("spectral_craddock26", n_clusters=100)

clu.plot_voxels_3d()

score = clu.report_silhouette()
```

---

# Architecture

## Layer 1 – Data Model

Handles input validation and standardized internal data representation.

Main functions:

* `make_fmri_tensor`
* `make_roi_table`
* `validate_data`
* `split_train_test`

Main class:

* `FMRIData`

---

## Layer 2 – Transform

Preprocessing and dimensionality reduction.

Implemented methods:

* Gaussian smoothing
* B-spline basis expansion
* Fourier basis expansion (planned)
* fPCA
* fICA (planned)
* Standardization

Main class:

* `PreprocessedData`

---

## Layer 3 – Similarity Computation

Generates similarity or distance matrices.

Supported metrics:

* Pearson correlation
* Partial correlation
* Cosine similarity
* L2 distance
* SRVF distance (planned)
* Affinity transformation (RBF, scaling)

Main class:

* `SimilarityResult`

---

## Layer 4 – Clustering

Graph-based and feature-based clustering.

Implemented methods:

* Spectral clustering (Craddock 26-neighbor constraint)
* Spectral clustering (unconstrained)
* K-means (planned)
* Hierarchical clustering
* Cluster postprocessing

Main class:

* `ClusteringResult`

---

## Layer 5 – Visualization & Reporting

Visualization:

* Similarity heatmaps
* Dendrogram
* 3D voxel scatter
* Cluster-colored voxel plots
* Figure saving utilities

Metrics:

* Silhouette score
* DICE coefficient

---

# Project Structure

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

* `data/` contains layer state classes.
* Other modules contain computational engines.
* `__init__.py` files expose a clean public API.

---

# Development Workflow

Recommended workflow:

1. Implement computational functions inside appropriate modules.
2. Connect them through layer class methods.
3. Commit frequently with meaningful messages.
4. Use feature branches for new algorithms.

---

# Roadmap

* Fourier basis expansion
* fICA implementation
* SRVF optimization
* GPU acceleration (future)
* Gradio GUI integration

---

# License

(To be specified)

---

# Author

Hyunseok Yoon

```

