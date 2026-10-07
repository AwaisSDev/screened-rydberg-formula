"""Reproduce the paper's worked examples and reference implementation (plain asserts; pytest optional).

    py -3.11 tests/test_paper_examples.py
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ionization import ie, ionization_energy  # noqa: E402


def close(a, b, tol):
    assert abs(a - b) <= tol, (a, b)


def test_oxygen_worked_example():          # paper section 4.7
    d = ie("O")
    close(d["sigma1"], 5.23348, 5e-6); close(d["T"], 3.3496, 5e-5); close(d["h"], 1.76652, 5e-6)
    close(d["D"], 0.71285, 5e-6); close(d["Zeff"], 2.05367, 5e-6); close(d["rydberg_eV"], 14.3457, 5e-5)
    close(d["hund_eV"], -0.8365, 5e-5); close(d["IE_eV"], 13.479, 5e-4)
    assert d["nu"] == {"same": 3, "in": 4, "core": 0, "df": 0, "out": 0}


def test_nitrogen_worked_example():        # paper section 4.7
    d = ie("N")
    close(d["sigma1"], 4.37867, 5e-5); close(d["Zeff"], 1.98350, 5e-5); close(d["hund_eV"], 0.8365, 5e-5)
    close(d["IE_eV"], 14.199, 5e-4)


def test_mg2plus():                        # paper section 4.7
    close(ionization_energy(12, 10), 78.63, 5e-3)


def test_grouping_rule_examples():         # paper section 2.2 / Table 6
    assert ie("Na")["nu"] == {"same": 0, "in": 8, "core": 2, "df": 0, "out": 0}
    assert ie("Fe")["nu"] == {"same": 1, "in": 14, "core": 10, "df": 0, "out": 0}   # 3d6 is in 'in', not 'df'
    assert ie("Pb")["nu"] == {"same": 1, "in": 20, "core": 60, "df": 0, "out": 0}
    close(ie("Pb")["sigma1"], 63.724, 5e-4)


def test_hydrogen_like():
    close(ionization_energy(1, 1), 13.598434, 1e-4)


def test_warnings():
    assert any("f-electron" in w for w in ie("Er", 66)["warnings"])
    assert any("Z > 110" in w for w in ie(118)["warnings"])
    assert any("heavy neutral" in w for w in ie("Pb")["warnings"])


def test_reference_implementation():       # tools/verify_from_inputs.py: paper equation + inputs table only
    r = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "verify_from_inputs.py")], capture_output=True,
                       text=True)
    assert r.returncode == 0, r.stdout + r.stderr


if __name__ == "__main__":
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            try:
                fn(); print("PASS", name)
            except AssertionError as e:
                fails += 1; print("FAIL", name, e)
    sys.exit(1 if fails else 0)
