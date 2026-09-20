# Future ML roadmap — planning only

No ML is implemented in this model case.

0. Current: establish X/M/Y schema, source provenance, independent structural views and explicit missingness.
1. Systematically collect matched in vivo datasets, prioritizing Tier-1 multi-drug studies.
2. Expand independent drug n; retain raw Clinical PT-level outcomes, provisional C1 grouping and common coverage rules. Do not replace missing PTs with a score.
3. Expand receptor representation (e.g. alpha3/alpha5) only after sufficient independent drugs and reproducible context exist; preserve receptor identity.
4. Compare prespecified predictive baselines: A Vina-only → Y; B PLIF-only → Y; C PLIF+Vina → Y; D PLIF+Vina+available M → Y. Decide outcomes, drug-level held-out evaluation, simple baseline, uncertainty and leakage controls before training. This future plan does not imply any current superiority or causal chain. Keep all poses/seeds of a drug in one split; fit preprocessing only on training drugs. Account for scaffold/study dependencies and reserve external drugs.
5. Investigate sparse-M approaches: multi-task, multi-view, missing-label/partially observed intermediate-phenotype methods and latent representations. Do not select algorithms before coverage and validation requirements are known; explicitly evaluate missingness assumptions.
6. External drug validation under a frozen protocol and endpoint definitions; report failures and calibration, not only successful correlations.
7. Only if validated, explore uncertainty/information-gain-based prioritization of new in vivo characterization. Active-learning or experimental-design criteria must be evaluated against cost, context compatibility and prospective validation.

n=4 is unsuitable for supervised model selection. Clinical reporting signals remain observational logROR, not incidence or risk ratios. An intermediate measured phenotype does not establish mediation. Sample size sufficiency must be evaluated for the eventual target/model and independent validation plan, not inferred from the 48 feature count.
