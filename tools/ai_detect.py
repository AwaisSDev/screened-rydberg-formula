"""Run the avoid-ai-writing detector (conorbronsdon/avoid-ai-writing) over a Markdown file, one
section at a time, because the bundled detector refuses documents longer than its size limit.

The detector is invoked as a subprocess so this project keeps no JavaScript dependency:
set AAW_REPO to a clone of the skill repository, e.g.
    git clone --depth 1 https://github.com/conorbronsdon/avoid-ai-writing %TEMP%\\aaw-repo

    py -3.11 tools/ai_detect.py docs/paper_draft.md            # per-section issues
    py -3.11 tools/ai_detect.py docs/paper_draft.md --summary  # one line per section
"""
import io
import json
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AAW = os.environ.get("AAW_REPO", os.path.join(tempfile.gettempdir(), "aaw-repo"))
CLI = os.path.join(AAW, "bin", "avoid-ai-writing.js")


def sections(path):
    """Split on '## ' headings; return (name, text) pairs."""
    out, name, buf = [], "(preamble)", []
    for line in io.open(path, encoding="utf-8").read().splitlines():
        if line.startswith("## "):
            out.append((name, "\n".join(buf)))
            name, buf = line[3:].strip(), []
        else:
            buf.append(line)
    out.append((name, "\n".join(buf)))
    return [(n, t) for n, t in out if len(re.findall(r"\w+", t)) > 40]


def analyse(text):
    fd, tmp = tempfile.mkstemp(suffix=".md")
    os.close(fd)
    io.open(tmp, "w", encoding="utf-8").write(text)
    try:
        p = subprocess.run(["node", CLI, "--context", "technical", tmp],
                           capture_output=True, text=True, encoding="utf-8", timeout=120)
        if not p.stdout.strip():
            return {"error": (p.stderr or "").strip()[:200]}
        return json.loads(p.stdout)
    finally:
        os.remove(tmp)


def main():
    path = sys.argv[1]
    summary = "--summary" in sys.argv
    if not os.path.exists(CLI):
        print("detector not found at", CLI, "- set AAW_REPO")
        return 1
    rows, total = [], 0
    for name, text in sections(path):
        r = analyse(text)
        n = len(r.get("issues", []))
        total += n
        rows.append((name, r.get("score"), r.get("label"), n, r))
    print("=" * 78)
    print("DETECTOR: %s" % os.path.relpath(path, ROOT))
    print("=" * 78)
    for name, score, label, n, r in rows:
        print("%-46s score=%-4s %-22s issues=%d" % (name[:46], score, str(label)[:22], n))
        if not summary:
            for it in r.get("issues", []):
                loc = it.get("line") or it.get("index") or ""
                msg = it.get("message") or it.get("description") or it.get("type") or ""
                txt = (it.get("text") or it.get("match") or "")[:70]
                print("      L%-4s %-28s %s | %s" % (loc, str(it.get("pattern") or it.get("category") or "")[:28],
                                                     msg[:60], txt))
    print("TOTAL ISSUES: %d over %d sections" % (total, len(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
