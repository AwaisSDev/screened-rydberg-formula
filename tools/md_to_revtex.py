r"""docs/paper_pra.md -> docs/pra/paper_pra.tex (REVTeX 4.2, Physical Review A, single-column 'preprint' layout).

    py -3.11 tools/to_pra.py && py -3.11 tools/md_to_revtex.py
Compile with pdflatex (e.g. upload the folder docs/pra/ to Overleaf; figures are copied next to the .tex).
"""
import io
import os
import re
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "docs", "paper_pra.md")
OUTDIR = os.path.join(ROOT, "docs", "pra")

UNI = {"§": r"\S", "²": r"\ensuremath{^{2}}", "³": r"\ensuremath{^{3}}", "·": r"\ensuremath{\cdot}",
       "¹": r"\ensuremath{^{1}}", "½": r"\ensuremath{\tfrac{1}{2}}", "×": r"\ensuremath{\times}",
       "Δ": r"\ensuremath{\Delta}", "Σ": r"\ensuremath{\Sigma}", "α": r"\ensuremath{\alpha}",
       "β": r"\ensuremath{\beta}", "δ": r"\ensuremath{\delta}", "ζ": r"\ensuremath{\zeta}", "κ": r"\ensuremath{\kappa}",
       "μ": r"\ensuremath{\mu}", "ν": r"\ensuremath{\nu}", "ρ": r"\ensuremath{\rho}", "σ": r"\ensuremath{\sigma}",
       "τ": r"\ensuremath{\tau}", "†": r"\ensuremath{^\dagger}", "…": r"\ldots{}", "′": r"\ensuremath{'}",
       "⁰": r"\ensuremath{^{0}}", "⁴": r"\ensuremath{^{4}}", "⁵": r"\ensuremath{^{5}}", "⁶": r"\ensuremath{^{6}}",
       "⁹": r"\ensuremath{^{9}}", "⁺": r"\ensuremath{^{+}}", "⁻": r"\ensuremath{^{-}}", "₀": r"\ensuremath{_{0}}",
       "₁": r"\ensuremath{_{1}}", "₂": r"\ensuremath{_{2}}", "→": r"\ensuremath{\to}", "−": r"\ensuremath{-}",
       "√": r"\ensuremath{\surd}", "∞": r"\ensuremath{\infty}", "≈": r"\ensuremath{\approx}",
       "≠": r"\ensuremath{\neq}", "≡": r"\ensuremath{\equiv}", "≤": r"\ensuremath{\leq}", "≥": r"\ensuremath{\geq}",
       "≳": r"\ensuremath{\gtrsim}", "⟨": r"\ensuremath{\langle}", "⟩": r"\ensuremath{\rangle}",
       "“": "``", "”": "''", "’": "'"}
PRE = r"""% single-column 'preprint' layout for submission; change preprint -> reprint for the two-column look
\documentclass[aps,pra,preprint,notitlepage,superscriptaddress,amsmath,amssymb]{revtex4-2}
\makeatletter
\@ifundefined{@titlepagefalse}{}{\@titlepagefalse}   % title, abstract and text on the same page
\makeatother
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{graphicx}
\usepackage{xurl}   % lets long URLs break anywhere
\emergencystretch=3em   % avoids lines running into the margin
\usepackage{newunicodechar}
\usepackage{adjustbox}   % shrinks over-wide tables to the text width
\usepackage{lineno}   % line numbers, as in the APS sample manuscript
%s
\begin{document}
\linenumbers
"""


def esc(t):
    t = re.sub(r"\\([*|_\[\]])", r"\1", t)          # Markdown escapes -> literal characters
    t = t.replace("\\", "\x00")
    for a, b in (("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"), ("}", r"\}"),
                 ("$", r"\$"), ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}"), ("|", r"\textbar{}")):
        t = t.replace(a, b)
    return t.replace("\x00", r"\textbackslash{}")


SUB = re.compile(r"(?<![\\A-Za-z_])(IE|[A-Za-zΑ-Ωα-ωΔΣ])_(\{[^}]*\}|[A-Za-z0-9′]+)(?![A-Za-z0-9_])")


def sub_tex(x, sub):
    """Z_eff -> $Z_{\mathrm{eff}}$, F_{n,j} -> $F_{n,j}$, r_c -> $r_{c}$."""
    if sub.startswith("{"):
        body = sub[1:-1]
    elif len(sub) > 1 and sub.isalpha():
        body = r"\mathrm{" + sub + "}"
    else:
        body = sub
    if x == "IE":
        x = r"\mathrm{IE}"
    return "$" + x + "_{" + body + "}$"


def inline(t, cite=True):
    parts = re.split(r"(`[^`]+`|https?://[^\s)]+)", t)
    out = []
    for p in parts:
        if p.startswith("`"):
            out.append(r"\texttt{" + esc(p[1:-1]) + "}")
        elif p.startswith("http"):
            tail = ""
            while p and p[-1] in ".,;:":
                tail, p = p[-1] + tail, p[:-1]
            out.append(r"\url{" + p + "}" + tail)
        else:
            subs = []

            def keep(m):
                if m.group(0) in ("d_near", "f_near"):      # screening-class names, not subscripts
                    return m.group(0)
                subs.append(sub_tex(m.group(1), m.group(2)))
                return "\x01%d\x02" % (len(subs) - 1)
            q = esc(SUB.sub(keep, p))
            q = re.sub("\x01(\\d+)\x02", lambda m: subs[int(m.group(1))], q)
            q = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", q)
            q = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", q)
            if cite:
                def ci(m):
                    nums = []
                    for a, b in re.findall(r"\[(\d+)\](?:–\[(\d+)\])?", m.group(0)):
                        nums += list(range(int(a), int(b) + 1)) if b else [int(a)]
                    return r"~\cite{" + ",".join(f"r{n}" for n in nums) + "}"
                q = re.sub(r" ?\[(?:[1-9]|[12]\d|3[0-4])\](?:(?:, |–)\[(?:[1-9]|[12]\d|3[0-4])\])*", ci, q)
            out.append(q)
    return "".join(out)


def table_tex(lines, caption=None):
    rows = [[c.strip() for c in re.split(r"(?<!\\)\|", l.strip()[1:-1])] for l in lines
            if not re.match(r"^\|[-| :]+\|$", l.strip())]
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]
    cells = [[inline(c, cite=False) for c in r] for r in rows]
    parts = [list(range(ncol))]
    if ncol > 9:    # very wide table: two blocks, each repeating the identifying first column(s)
        k = 2 if max(len(r[1]) for r in rows[1:]) <= 14 else 1
        rest = list(range(k, ncol))
        half = (len(rest) + 1) // 2
        parts = [list(range(k)) + rest[:half], list(range(k)) + rest[half:]]
    numeric = re.compile(r"^[\s\d.,+\-−–±%()/·†*\[\]a-z⁻⁰¹²³⁴⁵⁶⁷⁸⁹]*$")

    def col_align(c):   # numbers centred, text left
        vals = [r[c] for r in rows[1:] if r[c]]
        ok = lambda v: len(v) < 26 and (re.search(r"\d", v) or v.strip() in ("–", "-", "none")) \
            and numeric.match(v.replace("**", "")) or v.strip() == "none"
        return "c" if c > 0 and vals and all(ok(v) for v in vals) else "l"
    boxes = []
    for cols in parts:
        body = [" & ".join(r[c] for c in cols) + r" \\" for r in cells]
        tab = (r"\begin{tabular}{" + "".join(col_align(c) for c in cols) + "}\n\\hline\\hline\n" + body[0] + "\n\\hline\n"
               + "\n".join(body[1:]) + "\n\\hline\\hline\n" + r"\end{tabular}")
        # APS-style double rules; adjustbox shrinks a table only if it is wider than the text
        boxes.append(r"\begin{adjustbox}{max width=\linewidth}" + "\n" + tab + "\n" + r"\end{adjustbox}")
    box = "\n\\medskip\n".join(boxes)
    wide = ncol > 4 or max(len(" ".join(r)) for r in rows) > 70
    # single line spacing and a little column padding inside every table float; split tables are set smaller
    fmt = (r"\scriptsize" if len(parts) > 1 else r"\footnotesize") + \
        r"\linespread{1}\selectfont\setlength{\tabcolsep}{5pt}\renewcommand{\arraystretch}{1.1}"
    if caption is None:
        return "\\begin{table}[h]\n" + fmt + "\n" + box + "\n\\end{table}\n"
    env = "table*" if wide else "table"
    return (f"\\begin{{{env}}}[htbp]\n\\caption{{{inline(caption, cite=False)}}}\n{fmt}\n{box}\n"
            f"\\end{{{env}}}\n")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    s = io.open(SRC, encoding="utf-8").read()
    lines = s.split("\n")
    title = lines[0].lstrip("# ").strip()
    author = re.search(r"^\*\*(.+?)\*\*$", s, re.M).group(1)
    email = re.search(r"[\w.]+@[\w.]+\.\w+", s).group(0)
    i_ab = s.index("---") + 3
    i_ab_end = s.index("---", i_ab)
    abstract = " ".join(s[i_ab:i_ab_end].split())
    body_md = s[i_ab_end + 3:]
    i_ref = body_md.index("## REFERENCES")
    i_ref_end = body_md.index("---", i_ref)
    refs = re.findall(r"^\[(\d+)\] (.+)$", body_md[i_ref:i_ref_end], re.M)
    body_md = body_md[:i_ref] + body_md[i_ref_end + 3:]

    out = []
    L = body_md.split("\n")
    i, para, pending_caption, list_stack = 0, [], None, []

    def flush():
        nonlocal para
        if para:
            out.append(inline(" ".join(x.strip() for x in para)) + "\n")
            para = []

    def close_lists(level=0):
        while len(list_stack) > level:
            out.append("\\end{%s}" % list_stack.pop()[0])

    while i < len(L):
        ln = L[i]
        st = ln.strip()
        m_li = re.match(r"^(\s*)([-*]|\d+\.) (.*)$", ln)
        if not st or st == "---":
            flush()
            if not st and i + 1 < len(L) and re.match(r"^(\s*)([-*]|\d+\.) ", L[i + 1] or ""):
                i += 1
                continue
            close_lists()
            i += 1
            continue
        if st.startswith("$$"):
            flush(); close_lists()
            eq = [st]
            while not (eq[-1].endswith("$$") and (len(eq) > 1 or len(st) > 2 and st.endswith("$$") and st != "$$")):
                i += 1
                eq.append(L[i].strip())
            tex = "\n".join(eq).strip("$").strip()
            tex = re.sub(r"\s*\\tag\{[^}]*\}", "", tex)
            out.append("\\begin{equation}\n" + tex + "\n\\end{equation}")
            i += 1
            continue
        h = re.match(r"^(#{2,4}) (.+)$", st)
        if h:
            flush(); close_lists()
            lev, txt = len(h.group(1)), h.group(2)
            if lev == 2:
                if txt.startswith("ACKNOWLEDGMENTS"):
                    j = i + 1
                    ack = []
                    while j < len(L) and not L[j].startswith("## "):
                        if L[j].strip() and L[j].strip() != "---":
                            ack.append(L[j].strip())
                        j += 1
                    out.append("\\begin{acknowledgments}\n" + inline(" ".join(ack)) + "\n\\end{acknowledgments}")
                    i = j
                    continue
                if txt.startswith("APPENDIX"):
                    out.append("\\appendix\n\\section{" + inline(txt.split(":", 1)[1].strip().capitalize()) + "}")
                elif re.match(r"^[IVX]+\. ", txt):
                    out.append("\\section{" + inline(txt.split(". ", 1)[1].capitalize()) + "}")
                else:
                    out.append("\\section*{" + inline(txt.capitalize()) + "}")
            elif lev == 3:
                out.append("\\subsection{" + inline(re.sub(r"^[A-Z]\. ", "", txt)) + "}")
            else:
                out.append("\\subsubsection{" + inline(re.sub(r"^\d+\. ", "", txt)) + "}")
            i += 1
            continue
        img = re.match(r"^!\[[^\]]*\]\(([^)]+)\)$", st)
        if img:
            flush(); close_lists()
            src = os.path.normpath(os.path.join(ROOT, "docs", img.group(1)))
            shutil.copy(src, OUTDIR)
            j = i + 1
            while j < len(L) and not L[j].strip():
                j += 1
            cap = []
            if L[j].startswith("**FIG."):
                while j < len(L) and L[j].strip():
                    cap.append(L[j].strip()); j += 1
            capt = re.sub(r"^\*\*FIG\. \d+\.\*\* ?", "", " ".join(cap))
            out.append("\\begin{figure}[htbp]\n\\includegraphics[width=0.75\\columnwidth]{" + os.path.basename(src)
                       + "}\n\\caption{" + inline(capt) + "}\n\\end{figure}")
            i = j
            continue
        if st.startswith("**TABLE"):
            flush(); close_lists()
            cap = []
            while i < len(L) and L[i].strip():
                cap.append(L[i].strip()); i += 1
            pending_caption = re.sub(r"^\*\*TABLE [IVX]+\.\*\* ?", "", " ".join(cap))
            continue
        if st.startswith("|"):
            flush(); close_lists()
            tl = []
            while i < len(L) and L[i].strip().startswith("|"):
                tl.append(L[i]); i += 1
            out.append(table_tex(tl, pending_caption))
            pending_caption = None
            continue
        if m_li:
            flush()
            ind = len(m_li.group(1)) // 2
            kind = "enumerate" if m_li.group(2)[0].isdigit() else "itemize"
            while len(list_stack) > ind + 1:
                out.append("\\end{%s}" % list_stack.pop()[0])
            if len(list_stack) < ind + 1:
                out.append("\\begin{%s}" % kind); list_stack.append((kind, ind))
            item = [m_li.group(3)]
            i += 1
            while i < len(L) and L[i].strip() and not re.match(r"^(\s*)([-*]|\d+\.) ", L[i]) \
                    and L[i].startswith(" "):
                item.append(L[i].strip()); i += 1
            out.append("\\item " + inline(" ".join(item)))
            continue
        para.append(ln)
        i += 1
    flush(); close_lists()

    used = sorted({ch for ch in "".join(out) + title + abstract if ord(ch) > 127 and ch in UNI})
    nuc = "\n".join(f"\\newunicodechar{{{ch}}}{{{UNI[ch]}}}" for ch in used)
    bib = "\\begin{thebibliography}{%d}\n" % len(refs) + "\n".join(
        f"\\bibitem{{r{n}}} " + re.sub(r"QQCITE(\d+)QQ", lambda m: "Ref.~\\cite{r" + m.group(1) + "}",
                                       inline(re.sub(r"Ref\. (\d+)", r"QQCITE\1QQ", t), cite=False))
        for n, t in refs) + "\n\\end{thebibliography}"
    tex = (PRE.replace("%s", nuc) + f"\\title{{{inline(title)}}}\n\\author{{{author}}}\n\\email{{{email}}}\n"
           f"\\affiliation{{Independent researcher, Multan, Pakistan}}\n\\date{{\\today}}\n\\begin{{abstract}}\n{inline(abstract)}\n"
           "\\end{abstract}\n\\maketitle\n\n" + "\n\n".join(out) + "\n\n" + bib + "\n\\end{document}\n")
    path = os.path.join(OUTDIR, "paper_pra.tex")
    io.open(path, "w", encoding="utf-8", newline="\n").write(tex)
    leftover = sorted({ch for ch in tex if ord(ch) > 127 and ch not in UNI and ch not in "ÜáèéíöüĕńŠ–—"})
    print("wrote", path, "| refs", len(refs), "| unmapped non-ASCII:", leftover)


if __name__ == "__main__":
    main()
