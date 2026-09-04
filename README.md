# EEG Source Localization: Inverse Problem Solving via Regularization & Truncated SVD

![MATLAB](https://img.shields.io/badge/Language-MATLAB-orange.svg)
![Field](https://img.shields.io/badge/Domain-Biomedical%20Signal%20Processing%20%7C%20Inverse%20Problems-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

An academic computational project developed for the *Numerical Methods for Data Mining* course. This project models the electroencephalography (EEG) forward problem within a simplified spherical volume conductor and implements numerical regularization techniques to solve the severely ill-posed, underdetermined inverse problem of neural source localization.

---

## 📌 Executive Summary

Electroencephalography (EEG) source localization aims to estimate the non-invasive spatial distribution and temporal dynamics of intracranial neural current dipoles from scalp electric potential measurements. 

Because the number of internal active brain sources ($N = 32$) far exceeds the number of recorded scalp electrode channels ($M = 5$), the forward operator yields a severely underdetermined and ill-conditioned linear system:

$$V(t) = G x(t) + W(t), \quad M \ll N$$

This repository demonstrates:
1. **Geometric forward modeling**: Construction of a 3D spherical head volume conductor, random dipole distribution, and synthetic electrode placement.
2. **Dynamic time-series simulation**: Multichannel autoregressive (AR) signal generation modeling realistic cross-channel electrophysiological dependencies and additive Gaussian noise.
3. **Analytic minimum-norm estimation (MNE)**: Application of the right pseudo-inverse $x = G^T (GG^T)^{-1} V$.
4. **Truncated SVD (TSVD)**: Regularization via Moore-Penrose pseudo-inversion with singular-value thresholding ($	au \in [10^{-1}, 10^{-4}]$) to analyze stability, noise suppression, and rank truncation tradeoffs.

---

## 🧠 Mathematical & Numerical Formulation

### 1. Forward Modeling & Lead Field Matrix ($G$)
Under quasi-static approximations of Maxwell's equations, the forward mapping from source dipoles to surface potentials is linear:

$$V = G \cdot x$$

where:
* $V \in \mathbb{R}^{M \times T}$ represents the recorded potential at $M = 5$ electrodes over $T = 100$ time instances.
* $x \in \mathbb{R}^{N \times T}$ represents current dipole source activations for $N = 32$ distributed neural sources.
* $G \in \mathbb{R}^{M \times N}$ is the **Lead Field (Gain) Matrix**, computed here via an inverse-squared Euclidean distance metric representing volume conduction:

$$G_{ij} = \frac{c}{\|\mathbf{p}_{\text{sensor}, i} - \mathbf{p}_{\text{source}, j}\|^2}, \quad c = 1$$

### 2. Inverse Problem & Regularization Techniques

Since $M = 5 \ll N = 32$, the system has infinitely many solutions (null space $\mathcal{N}(G) \neq \{0\}$). Two approaches are implemented:

#### A. Closed-Form Minimum-Norm Estimate (MNE)
Selecting the unique solution that minimizes the Euclidean $\ell_2$-norm of source activity $\|x\|_2$:

$$\min_x \|x\|_2^2 \quad \text{subject to} \quad G x = V$$

Yielding the analytic right Moore-Penrose pseudo-inverse:

$$x_{\text{MNE}} = G^T (G G^T)^{-1} V$$

#### B. Truncated Singular Value Decomposition (TSVD)
Using the SVD of the gain matrix $G = U \Sigma V^T$, the Moore-Penrose pseudo-inverse $G^\dagger$ is computed by discarding singular values below a prescribed tolerance threshold $\tau$:

$$G^\dagger_\tau = \sum_{\sigma_i > \tau} \frac{1}{\sigma_i} v_i u_i^T$$

The code benchmarks the sensitivity of the reconstruction across four cutoff thresholds:
* $\tau = 10^{-1}$ (Over-regularized; excessive rank reduction, severe signal attenuation)
* $\tau = 10^{-2}$
* $\tau = 10^{-3}$ (Marginal threshold; noticeable reconstruction distortion)
* $\tau = 10^{-4}$ (Robust balance; preserves essential signal dynamics while discarding degenerate components)

Reconstruction discrepancy is quantitatively tracked via Frobenius / spectral norm metrics:

$$\epsilon_\tau = \|x_{\text{MNE}} - x_{\text{TSVD}}(\tau)\|$$

---

## 💻 Simulation Pipeline

```
  [ 3D Sphere Geometry ]        [ Autoregressive Dynamics ]
   R = 20, 32 Dipoles            5-Channel Coupled AR Model
           │                                 │
           ▼                                 ▼
   Lead Field Matrix G           Synthetic Potentials V(t)
     (5 sensors x 32 sources)       (5 channels x 100 time points)
           │                                 │
           └─────────────────┬───────────────┘
                             ▼
              [ Inverse Problem Solvers ]
            ┌─────────────────────────────┐
            │ 1. Minimum-Norm (Analytic)  │
            │ 2. Truncated SVD (4 Scales) │
            └──────────────┬──────────────┘
                           ▼
          [ Temporal & Error Analysis ]
          - Signal Tracking per Electrode
          - 2D Contour Error Residuals |V - G*x|
```

1. **Volume Conductor**: Sphere of radius $R = 20$ centered at the origin.
2. **Source Space**: 32 neural current sources randomly positioned in the interior ($r < 19$).
3. **Electrode Array**: 5 scalp electrodes distributed across the spherical surface.
4. **Time Series**: 100 discrete time instances simulated using coupled cross-channel delay terms and additive white Gaussian noise.
5. **Validation**: Direct comparison of ground truth scalp potentials $V(t)$ against re-projected estimates $\hat{V}(t) = G \hat{x}(t)$ across all 5 recording channels.

---

## 📊 Visualizations & Diagnostic Outputs

The script automatically generates comprehensive diagnostic plots:
* **Geometry Plots**: 3D spatial visualization of internal source locations within the translucent spherical volume conductor, alongside the surface electrode positions.
* **Electrode Time Tracking**: Multi-panel tiled layouts comparing:
  * Ground truth potentials $V_i(t)$
  * Analytic minimum-norm re-projection $G_i x(t)$
  * Stable regularized reconstruction with $\tau = 10^{-4}$
  * Degraded reconstruction under over-truncation with $\tau = 10^{-3}$
* **Residual Contour Maps**: Color-mapped error surfaces $|V - G \hat{x}|$ across space and time to assess spatiotemporal reconstruction bias.

---

## 🚀 Getting Started

### Prerequisites
* **MATLAB** (R2019b or newer recommended; requires `tiledlayout` support).
* No external toolboxes required (runs entirely on base MATLAB linear algebra routines).

### Running the Code
1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/eeg-source-localization-regularization.git
   cd eeg-source-localization-regularization
   ```
2. Open MATLAB, navigate to the cloned folder, and run:
   ```matlab
   run('eeg_source_localization.m')
   ```

---

## 🛠 Skills & Competencies Demonstrated

* **Numerical Linear Algebra**: Underdetermined linear systems, Moore-Penrose pseudo-inverses, Singular Value Decomposition (SVD), rank conditioning, and spectral matrix norms.
* **Regularization Theory**: Tikhonov / Minimum-norm formulation, Truncated SVD, noise sensitivity, and hyperparameter/threshold tuning.
* **Biomedical Modeling**: Biophysical volume conduction, forward/inverse EEG problem modeling, lead field matrix formulation.
* **MATLAB Computing**: Vectorized computation, 3D geometric visualization, dynamic simulation, and publication-ready multi-axis plotting.

---

## 👤 Author

**Chiara Giavalisco**  
* Master's Degree Coursework: *Metodi Numerici per il Data Mining*  
* [LinkedIn Profile]([https://linkedin.com](https://www.linkedin.com/in/chiara-giavalisco-28b1b9268/)) • [GitHub Profile]([https://github.com](https://github.com/chiaragiavalisco/chiaragiavalisco.github.io)) • [Email](mailto:chiara.giavalisco@gmail.com)

---

## 📄 License
This project is open-source and available under the [MIT License](LICENSE).
