"""Run the bundled avoid-ai-writing detector over a paper that is longer than
its MAX_WORDS (10000) ceiling, by splitting at section/paragraph boundaries.

The detector returns UNSCORED / "Text too long" for the whole document, so this
splits it into chunks under the limit, scores each one with the real detector
(node skills/ai-writing-detector/scripts/detect.js) and aggregates the results.

Detect-only: nothing is written except the temporary chunk files.

Usage:
  py -3.11 tools/ai_detect_chunked.py <input.txt|input.md> [--context technical]
                                     [--source-mode plain] [--limit 9000]
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

REPO = os.path.join(os.environ.get("TEMP", tempfile.gettempdir()), "aaw-repo")
DETECT = os.path.join(REPO, "skills", "ai-writing-detector", "scripts", "detect.js")


def read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def split_sections(text):
    """Split markdown at top-level headings, else at blank-line paragraphs."""
    lines = text.splitlines()
    heads = [
        i for i, ln in enumerate(lines)
        if re.match(r"^#{1,3}\s+\S", ln)
    ]
    if len(heads) >= 3:
        blocks, cur = [], []
        for i, ln in enumerate(lines):
            if i in heads and cur:
                blocks.append("\n".join(cur))
                cur = [ln]
            else:
                cur.append(ln)
        if cur:
            blocks.append("\n".join(cur))
    else:
        blocks = [b for b in re.split(r"\n\s*\n", text) if b.strip()]
    return [b for b in blocks if b.strip()]


def chunk(blocks, limit):
    """Greedy pack blocks into chunks of at most `limit` words."""
    chunks, cur, cur_words = [], [], 0
    for block in blocks:
        w = len(block.split())
        if w > limit:  # one oversized block: split on sentences
            if cur:
                chunks.append("\n\n".join(cur))
                cur, cur_words = [], 0
            words, piece = block.split(), []
            for sent in re.split(r"(?<=[.!?])\s+", block):
                sw = len(sent.split())
                if piece and len(" ".join(piece).split()) + sw > limit:
                    chunks.append(" ".join(piece))
                    piece = []
                piece.append(sent)
            if piece:
                chunks.append(" ".join(piece))
            continue
        if cur_words + w > limit and cur:
            chunks.append("\n\n".join(cur))
            cur, cur_words = [], 0
        cur.append(block)
        cur_words += w
    if cur:
        chunks.append("\n\n".join(cur))
    return chunks


def run_detector(text, context, source_mode, workdir, idx):
    path = os.path.join(workdir, "chunk_%02d.txt" % idx)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    cmd = ["node", DETECT, "--file", path,
           "--context", context, "--source-mode", source_mode]
    proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                          errors="replace")
    if proc.returncode not in (0, 2) or not proc.stdout.strip():
        return {"error": (proc.stderr or proc.stdout or "no output").strip()}
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        return {"error": "unparseable JSON: " + proc.stdout[:200]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("--context", default="technical")
    ap.add_argument("--source-mode", default="plain")
    ap.add_argument("--limit", type=int, default=9000)
    args = ap.parse_args()

    if not os.path.isfile(DETECT):
        sys.exit("detector not found at %s" % DETECT)

    text = read_text(args.input)
    total_words = len(text.split())
    chunks = chunk(split_sections(text), args.limit)
    workdir = tempfile.mkdtemp(prefix="aaw-chunks-")

    print("=" * 72)
    print("AI-WRITING DETECTOR (bundled v3.37.0) - chunked scan")
    print("file        : %s" % args.input)
    print("words       : %d   (detector ceiling is 10000, so it is split)" % total_words)
    print("chunks      : %d   context=%s  source-mode=%s" % (len(chunks), args.context, args.source_mode))
    print("=" * 72)

    results, tally = [], {}
    for i, ch in enumerate(chunks, 1):
        res = run_detector(ch, args.context, args.source_mode, workdir, i)
        wc = len(ch.split())
        if "error" in res:
            print("%2d. %5d words  ERROR: %s" % (i, wc, res["error"]))
            results.append({"words": wc, "error": res["error"]})
            continue
        issues = res.get("issues", []) or []
        for iss in issues:
            key = (iss.get("type") or iss.get("id") or "?",
                   iss.get("severity") or "?")
            tally[key] = tally.get(key, 0) + 1
        results.append({"words": wc, "score": res.get("score"),
                        "label": res.get("label"),
                        "classification": res.get("document_classification"),
                        "confidence": res.get("confidence_category"),
                        "issues": issues, "stats": res.get("stats", {}),
                        "too_long": res.get("tooLong", False)})
        print("%2d. %5d words  score=%-3s label=%-18s issues=%d%s"
              % (i, wc, res.get("score"), res.get("label"), len(issues),
                 "  [TOO LONG]" if res.get("tooLong") else ""))

    scored = [r for r in results if "score" in r and not r.get("too_long")]
    total_issues = sum(len(r.get("issues", [])) for r in scored)
    print("-" * 72)
    print("chunks scored        : %d" % len(scored))
    print("words scored         : %d" % sum(r["words"] for r in scored))
    print("total issues         : %d" % total_issues)
    print("issues per 1000 words: %.2f" % (1000.0 * total_issues / max(1, sum(r["words"] for r in scored))))
    if scored:
        avg = sum(r["score"] for r in scored) / len(scored)
        print("mean chunk score     : %.1f  (detector score is per-text, not additive)" % avg)
        print("worst chunk          : %s" % max(scored, key=lambda r: r["score"] or 0)["label"])

    print("-" * 72)
    print("ISSUES BY TYPE")
    if not tally:
        print("  (none)")
    for (typ, sev), n in sorted(tally.items(), key=lambda kv: -kv[1]):
        print("  %-34s %-9s %d" % (typ, sev, n))

    print("-" * 72)
    print("SEVERITY TOTALS")
    sev_tally = {}
    for (typ, sev), n in tally.items():
        sev_tally[sev] = sev_tally.get(sev, 0) + n
    for sev, n in sorted(sev_tally.items(), key=lambda kv: -kv[1]):
        print("  %-9s %d" % (sev, n))

    out = os.path.join(os.path.dirname(os.path.abspath(args.input)),
                       "aaw_chunked_report.json")
    try:
        with open(out, "w", encoding="utf-8") as fh:
            json.dump({"file": args.input, "words": total_words,
                       "chunks": results, "tally": {"%s|%s" % k: v for k, v in tally.items()}},
                      fh, indent=1)
        print("-" * 72)
        print("full JSON: %s" % out)
    except OSError as exc:
        print("could not write %s: %s" % (out, exc))


if __name__ == "__main__":
    main()
