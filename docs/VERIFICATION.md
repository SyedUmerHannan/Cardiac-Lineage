# Numerical verification

The main theorem is established by an analytical construction. The numerical tests check the implementation and provide reproducible joint-distribution witnesses.

| Check | Design |
|---|---|
| Reduced versus full formulation | 200 binary qualification masks, at the lower mean endpoint, midpoint and upper endpoint; 1,200 objective comparisons. Degenerate ranges include repeated mean values. |
| Constructive lift | 1,800 reconstructed laws checked for division and both daughter marginals, output mean, zero-output probability and nonnegative atom masses. |
| Two-stratum formula | 81 marginal pairs and five feasible means each; 405 comparisons with the closed-form frontier. |
| Boundary cases | Exact deadline equality, known failure atom, infeasible mean, invalid input and conditional-cohort validation. |
| Ranking reversal | Explicit paired-time distributions producing lower deadline output despite faster marginals. |
| Displayed frontier | 101 mean values with primal and dual diagnostics. |

The full and reduced numerical programmes use the same HiGHS solver with independently assembled constraints. They are formulation comparisons, not independent-solver validation. The two-stratum analytical formula and its explicit dual witness provide additional checks independent of numerical optimisation.

The constructive lift projects roundoff in status probabilities onto the admissible interval and normalises each row. Tests require nonnegative atom masses and allow numerical tolerance in marginal agreement. The LP uses floating-point arithmetic; diagnostic residuals are not interval-arithmetic certificates.

## Interpretation

The test inputs are synthetic. No biological parameters have been fitted. Six test methods cover the deadline model.

The exact theorem applies to finite-support population distributions. Midpoint binning of interval-censored timestamps can change deadline classification. Estimated marginals require a valid confidence set, and biological replication must respect donor or animal grouping. The observation schema specifies the required records.

## Related methods

The scalar sum problem has established fixed-marginal bounds (Makarov, DOI: 10.1137/1126086). Conditional Bernoulli bounds and linear-programming duality are standard ingredients. Transport methods already occur in cardiac and lineage analysis, including PILOT (DOI: 10.1038/s44320-023-00003-8) and moslin (DOI: 10.1186/s13059-024-03422-4).

The present formulation computes the complete shared-division-time mean–zero-output region and constructs compatible sister-pair laws. It does not reconstruct a preferred lineage map from molecular similarity. The biological use is checking compatibility of reported timing and endpoint measurements and determining which pairing information would identify the endpoints.
