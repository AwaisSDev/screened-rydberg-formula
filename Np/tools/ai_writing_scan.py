"""Scan a prose file for the AI-writing patterns catalogued by conorbronsdon/avoid-ai-writing.

The rule tables are parsed straight out of the skill's references/patterns.md so the scan uses the
published catalog rather than a hand-copied word list.

    py -3.11 tools/ai_writing_scan.py docs/paper_draft.md

Reports em dashes, Tier 1A/1B words, Tier 2 clusters, Tier 3 density, transition/template
phrases, bold density and paragraph-length uniformity. Counts only: no rewriting.
"""
import io
import os
import re
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "tools", "cache_patterns.md")
URL = "https://raw.githubusercontent.com/conorbronsdon/avoid-ai-writing/main/references/patterns.md"

EM_DASH = re.compile(r"\u2014|--")


def is_structure(line):
    """Markdown table separators ('|---|') and horizontal rules ('---') are layout, not
    punctuation. Counting their hyphen runs as em dashes reported 89 in a file that contains
    none, so they are filtered out of the dash count."""
    s = line.strip()
    return bool(s) and set(s) <= set("-:| \t")


def load_patterns():
    if os.path.exists(CACHE):
        return io.open(CACHE, encoding="utf-8").read()
    txt = urllib.request.urlopen(URL, timeout=30).read().decode("utf-8", "replace")
    io.open(CACHE, "w", encoding="utf-8").write(txt)
    return txt


def table_words(section, txt):
    """Pull the 'Replace' column out of the first markdown table after a heading."""
    i = txt.find(section)
    if i < 0:
        return []
    out, started = [], False
    for line in txt[i:].splitlines():
        if line.startswith("####") and started:
            break
        if line.startswith("|"):
            started = True
            cells = [c.strip() for c in line.strip("|").split("|")]
            if not cells or cells[0] in ("Replace", "Phrase") or set(cells[0]) <= set("-: "):
                continue
            cell = cells[0]
            for part in cell.split(" / "):
                part = re.sub(r"\(.*?\)", "", part).strip().lower()
                part = part.split()[0].strip("*`") if " " not in part else part.strip("*`")
                if part and len(part) > 2:
                    out.append(part)
        elif started and line.strip() == "":
            continue
    return out


def words(text):
    return re.findall(r"[A-Za-z][A-Za-z\u2019'-]*", text)


def scan(path):
    txt = io.open(path, encoding="utf-8").read()
    lines = txt.splitlines()
    pat = load_patterns()
    t1a = table_words("##### Tier 1A", pat)
    t1b = table_words("##### Tier 1B", pat)
    t2 = table_words("##### Tier 2", pat)
    t3 = table_words("##### Tier 3 -", pat) or table_words("##### Tier 3", pat)
    phrases = table_words("#### Tier 3 phrases", pat)

    trans = ["moreover", "furthermore", "additionally", "in today's", "in an era where",
             "it's worth noting", "notably", "in conclusion", "in summary", "to summarize",
             "when it comes to", "at the end of the day", "that said", "that being said",
             "here's what", "it is worth noting", "let's be clear", "to be honest"]

    print("=" * 72)
    print("FILE:", os.path.relpath(path, ROOT), " words:", len(words(txt)), " lines:", len(lines))
    print("=" * 72)

    dash_lines = "\n".join(ln for ln in lines if not is_structure(ln))
    n_dash = len(EM_DASH.findall(dash_lines))
    print("\n## Em dashes (target 0, hard max 1 per 1000 words): %d  [limit %d]"
          % (n_dash, max(1, len(words(txt)) // 1000)))
    for i, ln in enumerate(lines, 1):
        if not is_structure(ln) and EM_DASH.search(ln):
            print("  L%-4d %s" % (i, ln.strip()[:110]))

    def hits(label, terms):
        print("\n## %s" % label)
        total = 0
        for t in terms:
            rx = re.compile(r"\b" + re.escape(t) + r"\w*\b", re.I)
            found = [(i, rx.search(ln).group(0), ln.strip()[:90])
                     for i, ln in enumerate(lines, 1) if rx.search(ln)]
            if found:
                total += len(found)
                print("  %-20s %d  e.g. L%d %r" % (t, len(found), found[0][0], found[0][1]))
        print("  -- total %d" % total)
        return total

    hits("Tier 1A (AI frequency markers)", t1a)
    hits("Tier 1B (clarity edits)", t1b)
    hits("Tier 2 (flag in clusters)", t2)

    print("\n## Tier 3 (flag at >= max(3, 3%% of words) uses)")
    for t in t3:
        n = len(re.findall(r"\b" + re.escape(t) + r"\w*\b", txt, re.I))
        if n >= 3:
            print("  %-20s %d" % (t, n))

    print("\n## Tier 3 phrases (2+ uses, or 3+ distinct)")
    distinct = 0
    for p in phrases:
        n = len(re.findall(re.escape(p), txt, re.I))
        if n >= 2:
            print("  %-40s %d" % (p, n))
        if n:
            distinct += 1
    print("  distinct phrases present: %d (3+ = signal)" % distinct)

    print("\n## Transition / template phrases")
    for t in trans:
        n = len(re.findall(re.escape(t), txt, re.I))
        if n:
            print("  %-25s %d" % (t, n))

    body = [l for l in lines if l.strip() and not l.startswith(("#", "|", "![", "```", "$$", "-", "*", " "))]
    bold = sum(len(re.findall(r"\*\*[^*]+\*\*", l)) for l in body)
    print("\n## Bold in body prose: %d spans in %d paragraph lines" % (bold, len(body)))

    paras = [len(words(p)) for p in re.split(r"\n\s*\n", txt)
             if p.strip() and not p.lstrip().startswith(("#", "|", "!")
                                  ) and len(words(p)) > 20]
    if paras:
        print("## Paragraph lengths: n=%d min=%d max=%d mean=%.0f stdev=%.1f"
              % (len(paras), min(paras), max(paras), sum(paras) / len(paras),
                 (sum((x - sum(paras) / len(paras)) ** 2 for x in paras) / len(paras)) ** 0.5))
    return 0


if __name__ == "__main__":
    sys.exit(scan(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "docs", "paper_draft.md")))
