"""Check that every DOI printed in the paper resolves at doi.org to a live publisher page (HTTP status)."""
import re
import sys
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")
txt = open("docs/paper_pra.md", encoding="utf-8").read()
refs = txt[txt.index("## REFERENCES"):]
dois = sorted(set(re.findall(r"doi:\s*(10\.\d{4,9}/[^\s,;)]+?)(?:[.,;)]*(?:\s|$))", refs)))
dois = [d.rstrip(".") for d in dois]
print(len(dois), "DOIs")
bad = 0
for d in dois:
    url = "https://doi.org/" + d
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (ref-check)"}, method="HEAD")
        with urllib.request.urlopen(req, timeout=40) as r:
            print("200", d, "->", r.geturl()[:70])
    except urllib.error.HTTPError as e:
        # publishers (Elsevier, APS, Wiley) often answer 403/405 to scripts although the page exists
        print(e.code, d, "(publisher refused automated access; DOI itself resolved to" , e.url[:60] + ")")
        bad += e.code == 404
    except Exception as e:
        print("ERR", d, str(e)[:60])
        bad += 1
print("definitely broken (404/ERR):", bad)
