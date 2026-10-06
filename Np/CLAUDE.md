# Ionization-energy formula project

**Read `CONTINUE.md` first.** It has the goal, current status checklist, results so far, and exact resume steps. If the user says "continue", resume from the first unchecked item there.

- Run Python with `py -3.13` (the default `python` 3.14 has no packages). Don't create venvs.
  On the D:\ machine, which has no 3.13, use `py -3.11` with `export PYTHONPATH="D:/Chem-Research/Np/tools/numba_stub"`.
- `data/`, `common/` and `evaluate.py` are shared infrastructure; don't change their behaviour.
- Score any model with `py -3.13 evaluate.py results/<name>_predictions.csv` (columns Z,N,IE_pred_eV).
- Fitted models must report held-out results on the fixed splits: S1 test Z%5==0; S2 train Z<=54 / test Z>=55; S3 test N%6==0.
- Never SendMessage to an agent that is running inside a Workflow; it spawns a duplicate that clobbers files.
- Update `CONTINUE.md` (status checklist + timestamp) whenever a stage finishes.
