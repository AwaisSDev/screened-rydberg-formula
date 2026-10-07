"""Numeric-integrity check for a prose rewrite.

A style rewrite must not touch a single quantity. This compares the multiset of numeric tokens
(including superscripts, percentages, and decimal/coefficient forms) between the pre-rewrite backup
and the current file, and lists any token that appeared or disappeared.

    py -3.11 tools/ai_number_check.py docs/review/paper_before_ai_rewrite.md docs/paper_draft.md
"""
import io
import os
import re
import sys
from collections import Counter

TOKEN = re.compile(r"[0-9\u2070-\u2079\u00b2\u00b3\u00b9]+(?:[.,][0-9]+)*")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def tokens(path):
    txt = io.open(path, encoding="utf-8", newline="").read()
    return Counter(TOKEN.findall(txt)), len(re.findall(r"\w+", txt))


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:
        pass
    a, b = sys.argv[1], sys.argv[2]
    ta, wa = tokens(a)
    tb, wb = tokens(b)
    print("=" * 72)
    print("NUMERIC INTEGRITY")
    print("  before: %s  words=%d  numeric tokens=%d" % (os.path.relpath(a, ROOT), wa, sum(ta.values())))
    print("  after : %s  words=%d  numeric tokens=%d" % (os.path.relpath(b, ROOT), wb, sum(tb.values())))
    gone = ta - tb
    added = tb - ta
    if not gone and not added:
        print("  RESULT: every numeric token is preserved, none added.")
        return 0
    if gone:
        print("  LOST (%d distinct):" % len(gone))
        for k, v in sorted(gone.items()):
            print("    %-24s x%d" % (k, v))
    if added:
        print("  NEW (%d distinct):" % len(added))
        for k, v in sorted(added.items()):
            print("    %-24s x%d" % (k, v))
    return 1


if __name__ == "__main__":
    sys.exit(main())
