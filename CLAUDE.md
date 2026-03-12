# multiGGMr - Bayesian Inference of Multiple Gaussian Graphical Models

R package implementing Peterson, Stingo & Vannucci (2015, JASA) for joint inference
of multiple GGMs with shared structure via MRF priors. Reference MATLAB code:
https://github.com/cbpeterson/multiGGM

## Method Summary (Peterson et al. 2015)

**Goal**: Infer K undirected graphical models (precision matrices Ω_1,...,Ω_K) for K
related sample groups, borrowing strength across groups when appropriate.

### Model Components

1. **Likelihood**: X_k ~ N(μ_k, Ω_k^{-1}), k=1,...,K. Column-centered data → S_k = X_k'X_k.
2. **G-Wishart prior**: Ω_k | G_k ~ W_G(b, D). Conjugate. Default: b=3, D=I_p.
3. **MRF prior on graphs** (eq 3.2):
   P(g_ij | ν, Θ) ∝ exp(ν_ij * 1'g_ij + g_ij' Θ g_ij)
   where g_ij = (g_{1,ij},...,g_{K,ij})' is binary vector of edge inclusion across K groups.
   Normalizing constant C(Θ, ν_ij) = Σ_{g∈{0,1}^K} exp(ν_ij * 1'g + g'Θg). Requires 2^K terms.
4. **Spike-and-slab on Θ** (network similarity):
   θ_km | γ_km ~ (1-γ_km)δ_0 + γ_km * Gamma(α, β), α>1 (non-local).
   γ_km ~ Bernoulli(w). Paper defaults: α=2, β=5(rate), w=0.9.
5. **Edge-specific prior on ν**: q_ij = logit^{-1}(ν_ij) ~ Beta(a, b).
   Paper default: a=1, b=4 → prior P(edge) ≈ 0.20.

### MCMC Algorithm (Appendix A)

Per iteration t:
- **Step 1**: Update (G_k, Ω_k) for each k using Wang & Li (2012) exchange algorithm.
  For each edge (i,j): Bernoulli proposal from MRF conditional (eq A.3), then
  acceptance via ratio r2 of G-Wishart NO_{ij} pdfs (Step 2b), then
  posterior update of C given current graph (Step 2c), then
  BIPS clique update for full precision (Step 3).
- **Step 2**: Update (θ_km, γ_km) via MH between-model + within-model moves (eqs A.5-A.7).
  Between: toggle spike↔slab. Within: propose from Gamma(α_prop, β_prop).
- **Step 3**: Update ν_ij via MH with Beta(a_prop, b_prop) proposal (eqs A.8-A.10).

### Paper Hyperparameter Defaults (Section 5.1)
- G-Wishart: b=3, D=I_p
- Θ slab: Gamma(α=2, rate=5) → mean=0.4, P(θ≤1)=0.96
- γ prior: w=0.9 (strong belief networks related)
- ν prior: Beta(a=1, b=4) → prior edge prob ≈ 0.20
- ν proposal: Beta(a_prop=2, b_prop=4) → avg acceptance ~39%
- MCMC: 10K burnin, 20K post-burnin samples

## Package Architecture

### R files
- `R/multiggm_mcmc.R`: Main entry + single-chain MCMC engine (`.multiggm_mcmc_single`)
- `R/methods.R`: S3 print/summary methods
- `R/pcor.R`: Partial correlation utilities, PIP, ROC/AUC, confusion metrics
- `R/calc_mrf_logC.R`: R wrapper for C++ MRF normalizing constant
- `R/fix_matrix.R`: Matrix positive-definiteness helper

### C++ files (src/)
- `src/calc_mrf_logC.cpp`: MRF log-normalizing constant via exact enumeration over {0,1}^K
- `src/wangli.cpp`: Wang & Li (2012) G-Wishart samplers:
  - `GWishart_NOij_Gibbs_cpp`: Edge-specific Gibbs update
  - `GWishart_BIPS_maximumClique_cpp`: Block update via clique decomposition
  - `log_H_cpp`: Log acceptance ratio for edge proposal (eq A.3)
  - `log_GWishart_NOij_pdf_cpp`: G-Wishart pdf for exchange step
  - Helper: `wishrnd_cpp`, `log_iwishart_invA_const_cpp`, `log_J_cpp`

## Known Issues & Discrepancies

### Hyperparameter defaults don't match paper
- Code: alpha=2, beta=2 (Gamma scale=2 → mean=4). Paper: rate=5 → scale=0.2 → mean=0.4
- Code: w=0.1. Paper: w=0.9
- Code: a=1, b=1 (uniform q_ij). Paper: a=1, b=4
- Code: alpha_prop=1, beta_prop=1. Paper: not specified for Gamma proposal, but a_prop=2, b_prop=4 for Beta

### Performance bottleneck
- Main MCMC loop is R-level triple-nested: iter × p(p-1)/2 edges × K groups
- Each edge calls C++ functions individually → massive R↔C++ overhead
- The Theta/nu update loops also iterate over all edges in R calling calc_mrf_logC per edge
- **Fix**: Move entire per-iteration logic into C++, only return saved samples to R

### Code quality
- Extensive debug counters and commented-out code blocks
- message() on every iteration for nu acceptance floods output
- Several TODO comments indicating uncertainty about correctness
- Redundant bookkeeping variables

### Missing functionality
- No `simulate_data()` for generating synthetic test data
- No acceptance of raw data (only covariance matrices S_list)
- No progress reporting / progress bar
- No convergence diagnostics (R-hat, ESS)
- No vignette or worked example
- Several utility functions not exported (posterior_ci, pip_edges, etc.)
