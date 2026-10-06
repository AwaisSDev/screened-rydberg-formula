"""Track B (data-driven / pattern recognition) models for successive ionization energies.

Modules
-------
configs    : configuration helpers (removed-subshell logic, jj filling, Hund pair counts)
baselines  : Bohr, Slater (1-electron and total-energy), Clementi-Raimondi (both variants)
gshm       : the new Generalized Screened-Hydrogenic Model (predict(Z, N, shells=None))
fit_gshm   : fitting / held-out validation / sparse term selection / ablation (script)
explore    : pattern-recognition figures (script)
"""
