"""Read-only import of models/push_a/model.py under a unique module name (push_b also has a model.py).

    from push_a_import import PA      # PA.fit, PA.evaluate, PA.predict_many, PA.default_spec, ...
"""
import importlib.util
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_PATH = os.path.join(_ROOT, "models", "push_a", "model.py")

if "push_a_model" in sys.modules:
    PA = sys.modules["push_a_model"]
else:
    _spec = importlib.util.spec_from_file_location("push_a_model", _PATH)
    PA = importlib.util.module_from_spec(_spec)
    sys.modules["push_a_model"] = PA
    _spec.loader.exec_module(PA)
