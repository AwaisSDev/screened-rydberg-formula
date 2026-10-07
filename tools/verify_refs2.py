"""Second pass: details Crossref's search did not settle, and what each cited paper actually says (abstracts)."""
import json
import re
import sys
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "ref-check/1.0 (mailto:mawais9171@gmail.com)"}


def get(url, tries=3):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=40) as r:
                return json.load(r)
        except Exception as e:
            err = e
            time.sleep(2 + 3 * k)
    raise err


sys.stdout.reconfigure(encoding="utf-8")
print("=== [31] Yerokhin & Shabaev: article number")
m = get("https://api.crossref.org/works/10.1063/1.4927487")["message"]
print({k: m.get(k) for k in ("article-number", "volume", "issue", "page", "title")})

print("\n=== [8] Edlen: parent volume")
for doi in ("10.1007/978-3-662-35391-2", "10.1007/978-3-662-35391-2_2"):
    m = get("https://api.crossref.org/works/" + doi)["message"]
    print(doi, {k: m.get(k) for k in ("title", "container-title", "volume", "editor", "publisher", "issued", "page", "type")})

print("\n=== [25] Crilly: volume/article")
m = get("https://api.crossref.org/works/10.1016/j.hedp.2023.101053")["message"]
print({k: m.get(k) for k in ("volume", "article-number", "issued", "author")} if False else
      {k: m.get(k) for k in ("volume", "article-number", "issued")}, [a.get("family") for a in m.get("author", [])][:4])

print("\n=== [19] Di Rocco 1992: how do the citing papers print it?")
for doi in ("10.1590/s0103-97332005000100008", "10.1007/s13538-015-0397-9", "10.1016/j.hedp.2015.08.002"):
    m = get("https://api.crossref.org/works/" + doi)["message"]
    for ref in m.get("reference", []):
        txt = json.dumps(ref, ensure_ascii=False)
        if re.search(r"Di ?Rocco|Rocco", txt) and re.search(r"199[0-9]|1992", txt):
            print(doi, "->", txt[:300])

print("\n=== abstracts (Semantic Scholar) for what each paper is cited for")
DOIS = {4: "10.1016/0003-4916(59)90023-5", 5: "10.1016/0003-4916(64)90192-7", 6: "10.1103/RevModPhys.35.436",
        7: "10.1098/rspa.1958.0182", 9: "10.1088/0031-8949/47/3/007", 12: "10.1016/S0022-4073(97)00018-6",
        14: "10.1016/S0022-4073(97)00226-4", 15: "10.1016/S0022-4073(01)00142-X", 16: "10.1016/j.hedp.2011.04.006",
        17: "10.1088/0031-8949/29/5/005", 20: "10.1590/S0103-97332005000100008", 22: "10.1007/s13538-015-0397-9",
        23: "10.1103/PhysRevA.5.1137", 24: "10.1016/j.hedp.2005.07.001", 25: "10.1016/j.hedp.2023.101053",
        28: "10.1103/PhysRevA.55.191", 29: "10.1103/PhysRevA.47.3649", 30: "10.1016/j.adt.2003.11.005",
        31: "10.1063/1.4927487", 21: "10.1016/j.hedp.2015.08.002"}
for n, doi in DOIS.items():
    try:
        d = get("https://api.semanticscholar.org/graph/v1/paper/DOI:" + urllib.parse.quote(doi) + "?fields=title,abstract")
        ab = (d.get("abstract") or "").replace("\n", " ")
        print(f"[{n}] {d.get('title','')[:80]}\n     {ab[:420] if ab else '(no abstract in Semantic Scholar)'}")
    except Exception as e:
        print(f"[{n}] Semantic Scholar: {str(e)[:60]}")
    time.sleep(1.2)
