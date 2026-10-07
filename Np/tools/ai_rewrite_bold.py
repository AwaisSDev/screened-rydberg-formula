"""Strip emphasis bold from body prose in the paper draft.

The avoid-ai-writing catalog ("Bold overuse") asks for at most one bolded phrase per major section,
which for this manuscript means none in flowing prose. Table cells, table and figure caption labels
("**Table 1.**", "**Figure 2.**") and the abstract/index-terms labels are structural, so they keep
their bold. Display equations, code spans and reference entries are left untouched.

    py -3.11 tools/ai_rewrite_bold.py docs/paper_draft.md [--apply]

Without --apply it only reports what it would change.
"""
import io
import re
import sys

KEEP = re.compile(r"^\s*\*\*(Table|TABLE|Figure|FIG|Appendix|Abstract|Index Terms)[ .]")
BOLD = re.compile(r"\*\*(.+?)\*\*")


def strip_line(ln):
    """Remove ** emphasis from one prose line, keeping any table-cell wording intact."""
    if ln.lstrip().startswith(("|", "#", "![", "$$", ">", "```")) or KEEP.match(ln):
        return ln, []
    found = [m.group(1) for m in BOLD.finditer(ln)]
    if not found:
        return ln, []
    return BOLD.sub(r"\1", ln), found


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:
        pass
    path = sys.argv[1]
    apply_ = "--apply" in sys.argv
    # newline="" keeps the original CRLF endings: each line keeps its trailing "\r",
    # so joining on "\n" reproduces the file byte for byte apart from the edits.
    lines = io.open(path, encoding="utf-8", newline="").read().split("\n")
    out, changes, in_eq, in_code = [], 0, False, False
    for i, ln in enumerate(lines, 1):
        if ln.strip().startswith("$$"):
            in_eq = not in_eq
            out.append(ln)
            continue
        if ln.strip().startswith("```"):
            in_code = not in_code
            out.append(ln)
            continue
        if in_eq or in_code:
            out.append(ln)
            continue
        new, found = strip_line(ln)
        if found:
            changes += len(found)
            print("L%-4d -%s" % (i, " -".join(found)[:100]))
        out.append(new)
    print("bold spans removed from body prose: %d" % changes)
    if apply_:
        io.open(path, "w", encoding="utf-8", newline="\n").write("\n".join(out))
        print("wrote", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
