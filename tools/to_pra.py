"""Convert docs/paper_draft.md (IEEE-style draft) to Physical Review A conventions -> docs/paper_pra.md.

    py -3.11 tools/to_pra.py && py -3.14 tools/build_paper_pdf.py docs/paper_pra.md

APS / Physical Review conventions applied (APS Style and Notation Guide; PRA Information for Authors):
  * sections I., II., ... in capitals; subsections A., B.; sub-subsections 1., 2.; cross references "Sec. II B";
  * numbered display equations (1), (2), ...; appendix equations (A1); text references "Eq. (3)";
  * table captions "TABLE I." (Roman), figure captions "FIG. 1."; text "Table I", "Fig. 1";
  * references in order of citation, Physical Review style with titles: A. B. Author, Title, Journal **vol**, page (year);
  * no keyword list (PRA uses PhySH subject headings chosen at submission); abstract without label;
  * AI-assistance disclosure placed in the Acknowledgments (APS policy on appropriate use of AI tools).
"""
import io
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "docs", "paper_draft.md")
DST = os.path.join(ROOT, "docs", "paper_pra.md")
ROM = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII"]
ABC = "ABCDEFGHIJ"


PR_FIXED = {   # entries that are not journal articles: Physical Review style written out
    "B. Edlén": "B. Edlén, Atomic spectra, in *Handbuch der Physik*, Vol. 27, *Spectroscopy I*, edited by S. Flügge "
                "(Springer, Berlin, 1964), pp. 80–220, doi:10.1007/978-3-662-35391-2_2.",
    "H. Mayer": "H. Mayer, *Methods of Opacity Calculations*, Los Alamos Scientific Laboratory Report No. LA-647 (1947).",
    "R. D. Cowan": "R. D. Cowan, *The Theory of Atomic Structure and Spectra* (University of California Press, Berkeley, "
                   "1981).",
    "C. Froese Fischer, T. Brage": "C. Froese Fischer, T. Brage, and P. Jönsson, *Computational Atomic Structure: An MCHF "
                                   "Approach* (Institute of Physics Publishing, Bristol, 1997).",
    "A. Kramida": "A. Kramida, Yu. Ralchenko, J. Reader, and NIST ASD Team, NIST Atomic Spectra Database (version 5.12) "
                  "(National Institute of Standards and Technology, Gaithersburg, MD, 2024), https://physics.nist.gov/asd, "
                  "doi:10.18434/T4W30F, accessed on or before 5 October 2026.",
}


def pr_reference(s):
    """IEEE entry -> Physical Review style."""
    s = s.strip()
    for k, v in PR_FIXED.items():
        if s.startswith(k):
            m = re.search(r"as cited in \[(\d+)\]", s)
            return v.replace("[@P@]", m.group(1) if m else "")
    m = re.match(r"(?P<au>.+?), “(?P<ti>.+?),” \*(?P<jo>[^*]+)\*, vol\. (?P<vo>[^,]+),(?: no\. [^,]+,)? "
                 r"(?:pp\. |p\. )?(?P<pg>[A-Z]?\d+)(?:–[A-Z]?\d+)?, (?P<yr>\d{4})(?:, doi: (?P<doi>.+?))?\.$", s)
    if m:
        au = m["au"].replace(", and ", ", and ")
        out = f"{au}, {m['ti']}, {m['jo']} **{m['vo']}**, {m['pg']} ({m['yr']})"
        return out + (f", doi:{m['doi']}." if m["doi"] else ".")
    m = re.match(r"(?P<au>.+?), “(?P<ti>.+?),” \*(?P<jo>[^*]+)\*, vol\. (?P<vo>[^,]+), (?P<pg>\d+), (?P<yr>\d{4})\.$", s)
    if m:   # article-number style (J. Phys. Chem. Ref. Data 44, 033103)
        return f"{m['au']}, {m['ti']}, {m['jo']} **{m['vo']}**, {m['pg']} ({m['yr']})."
    m = re.match(r"(?P<au>.+?), “(?P<ti>.+?),” \*(?P<jo>[^*]+)\*, (?P<yr>\d{4}), doi: (?P<doi>.+?)\.$", s)
    if m:
        return f"{m['au']}, {m['ti']}, {m['jo']} ({m['yr']}), doi:{m['doi']}."
    return s   # books, reports, database entries: already in a readable form


def main():
    s = io.open(SRC, encoding="utf-8").read()
    # ---- front matter: abstract label off, keywords off
    s = s.replace("**Abstract—**", "", 1)
    s = re.sub(r"\n\*\*Index Terms—\*\*[^\n]*(\n[^\n]+)*?\n\n", "\n\n", s, count=1)
    # ---- sections
    sec, sub = {}, {}
    def h2(m):
        num, title = m.group(1), m.group(2)
        sec[num] = ROM[int(num)]
        return f"## {ROM[int(num)]}. {title.upper()}"
    s = re.sub(r"^## (\d+)\. (.+)$", h2, s, flags=re.M)
    def h3(m):
        a, b, title = m.group(1), m.group(2), m.group(3)
        sub[f"{a}.{b}"] = f"{ROM[int(a)]} {ABC[int(b) - 1]}"
        return f"### {ABC[int(b) - 1]}. {title}"
    s = re.sub(r"^### (\d+)\.(\d+) (.+)$", h3, s, flags=re.M)
    def h4(m):
        a, b, c, title = m.groups()
        sub[f"{a}.{b}.{c}"] = f"{ROM[int(a)]} {ABC[int(b) - 1]} {c}"
        return f"#### {c}. {title}"
    s = re.sub(r"^#### (\d+)\.(\d+)\.(\d+) (.+)$", h4, s, flags=re.M)
    s = re.sub(r"^## (Data and code availability|Acknowledgements|References|AI-assistance disclosure)$",
               lambda m: "## " + m.group(1).upper().replace("ACKNOWLEDGEMENTS", "ACKNOWLEDGMENTS"), s, flags=re.M)
    s = s.replace("## Appendix A. Fitted parameters (all-data fit)", "## APPENDIX: FITTED PARAMETERS (ALL-DATA FIT)")

    # ---- section cross references (skip references to other documents, e.g. `docs/x.md` §2.2)
    def secref(m):
        pre, a, b, c, d = m.group(1), m.group(2), m.group(3), m.group(4), m.group(5)
        if pre.rstrip().endswith(("md`", ".md", "md`)")):
            return m.group(0)
        key = a + (f".{b}" if b else "") + (f".{c}" if c else "")
        lab = sub.get(key) or sec.get(a)
        if d:   # range like 2.2–2.4
            return pre + f"Secs. {lab}–{sub.get(a + '.' + d, d)}"
        return pre + ("Sec. " + lab if lab else m.group(0)[len(pre):])
    s = re.sub(r"(.{0,6})§(\d)(?:\.(\d))?(?:\.(\d))?(?:–\d\.(\d))?", secref, s)

    # ---- equations: number displays, map 'eq. (2.2)' (the model equation) to its number
    n_eq, out, app = 0, [], False
    for block in re.split(r"(\$\$.*?\$\$)", s, flags=re.S):
        if block.startswith("$$"):
            if app:
                n_eq_lab = f"A{(n_eq := n_eq + 1) - main_count}"
            else:
                n_eq += 1
                n_eq_lab = str(n_eq)
            block = block[:-2].rstrip() + rf" \tag{{{n_eq_lab}}}" + "\n$$" if "\n" in block else \
                block[:-2] + rf" \tag{{{n_eq_lab}}}$$"
        elif "## APPENDIX" in block and not app:
            app, main_count = True, n_eq
        out.append(block)
    s = "".join(out)
    s = re.sub(r"\beq\. \(2\.2\)", "Eq. (3)", s)

    # ---- tables and figures
    tmap = {"1": "I", "2": "II", "3": "III", "4": "IV", "4a": "V", "5": "VI", "6": "VII", "A1": "VIII", "A2": "IX"}
    s = re.sub(r"\*\*Table (4a|A1|A2|\d)\.\*\*", lambda m: f"**TABLE {tmap[m.group(1)]}.**", s)
    s = re.sub(r"\bTables? (4a|A1|A2|\d)\b", lambda m: m.group(0).split()[0] + " " + tmap[m.group(1)], s)
    s = re.sub(r"Tables (\w+)–(\d)", lambda m: m.group(0), s)
    s = re.sub(r"\*\*Figure (\d)\.\*\*", r"**FIG. \1.**", s)
    s = re.sub(r"\bFigures? (\d)\b", r"Fig. \1", s)
    s = re.sub(r"!\[Fig\. (\d)\]", r"![Figure \1]", s)

    # ---- AI disclosure into the Acknowledgments
    m = re.search(r"## AI-ASSISTANCE DISCLOSURE\n\n(.*?)\n\n## ACKNOWLEDGMENTS\n\n(.*?)\n", s, flags=re.S)
    if m:
        ai = " ".join(m.group(1).split())
        s = s.replace(m.group(0), "## ACKNOWLEDGMENTS\n\n" + ai + "\n")

    # ---- references: Physical Review style
    i = s.index("## REFERENCES")
    j = s.index("---", i)
    refs = s[i:j]
    refs = re.sub(r"^\[(\d+)\] (.+)$", lambda m: f"[{m.group(1)}] " + pr_reference(m.group(2)), refs, flags=re.M)
    s = s[:i] + refs + s[j:]
    io.open(DST, "w", encoding="utf-8", newline="\n").write(s)
    print("wrote", DST, "| equations:", n_eq, "| sections:", sec)


if __name__ == "__main__":
    main()
