"""Check every article reference of the paper against Crossref (https://api.crossref.org).

    python tools/verify_refs.py          -> prints, per reference, what the paper says vs what Crossref returns

Nothing is assumed: a reference is only 'OK' when Crossref returns a record whose first-author family name, volume,
first page and year all agree with the paper. DOIs printed in the paper are fetched directly; for references
without a printed DOI the best Crossref hit for title + year is shown and compared the same way.
"""
import json
import re
import sys
import urllib.parse
import urllib.request

# n, first-author family, volume, first page, year, journal fragment, DOI printed in the paper (or None), title for search
REFS = [
    (1, "Slater", "36", "57", 1930, "Physical Review", None, "Atomic Shielding Constants"),
    (2, "Clementi", "38", "2686", 1963, "Chemical Physics", None, "Atomic Screening Constants from SCF Functions"),
    (3, "Clementi", "47", "1300", 1967, "Chemical Physics", None, "Atomic Screening Constants from SCF Functions. II. Atoms with 37 to 86 Electrons"),
    (4, "Layzer", "8", "271", 1959, "Annals of Physics", "10.1016/0003-4916(59)90023-5", None),
    (5, "Layzer", "29", "101", 1964, "Annals of Physics", "10.1016/0003-4916(64)90192-7", None),
    (6, "Scherr", "35", "436", 1963, "Reviews of Modern Physics", "10.1103/RevModPhys.35.436", None),
    (7, "Dalgarno", "247", "245", 1958, "Proceedings of the Royal Society", "10.1098/rspa.1958.0182", None),
    (8, "Edl", None, "80", 1964, "Handbuch", "10.1007/978-3-662-35391-2_2", None),
    (9, "Safronova", "47", "364", 1993, "Physica Scripta", "10.1088/0031-8949/47/3/007", None),
    (11, "More", "27", "345", 1982, "Quantitative Spectroscopy", "10.1016/0022-4073(82)90127-3", None),
    (12, "Faussurier", "58", "233", 1997, "Quantitative Spectroscopy", "10.1016/S0022-4073(97)00018-6", None),
    (13, "Faussurier", "4", "114", 2008, "High Energy Density Physics", None, "Equation of state of dense plasmas using a screened-hydrogenic model with l-splitting"),
    (14, "Martel", "60", "623", 1998, "Quantitative Spectroscopy", "10.1016/S0022-4073(97)00226-4", None),
    (15, "Rubiano", "72", "575", 2002, "Quantitative Spectroscopy", "10.1016/S0022-4073(01)00142-X", None),
    (16, "Mendoza", "7", "169", 2011, "High Energy Density Physics", "10.1016/j.hedp.2011.04.006", None),
    (17, "Kregar", "29", "438", 1984, "Physica Scripta", "10.1088/0031-8949/29/5/005", None),
    (18, "Kregar", "31", "246", 1985, "Physica Scripta", "10.1088/0031-8949/31/4/005", None),
    (20, "Pomarico", "35", "130", 2005, "Brazilian Journal of Physics", "10.1590/S0103-97332005000100008", None),
    (21, "Lanzini", "17", "240", 2015, "High Energy Density Physics", "10.1016/j.hedp.2015.08.002", None),
    (22, "Di Rocco", "46", "175", 2016, "Brazilian Journal of Physics", "10.1007/s13538-015-0397-9", None),
    (23, "Rozsnyai", "5", "1137", 1972, "Physical Review A", "10.1103/PhysRevA.5.1137", None),
    (24, "Chung", "1", "3", 2005, "High Energy Density Physics", None, "FLYCHK: Generalized population kinetics and spectral model for rapid spectroscopic analysis for all elements"),
    (25, "Crilly", None, None, 2023, "High Energy Density Physics", "10.1016/j.hedp.2023.101053", None),
    (26, "Koopmans", "1", "104", 1934, "Physica", None, "Über die Zuordnung von Wellenfunktionen und Eigenwerten zu den einzelnen Elektronen eines Atoms"),
    (27, "Kohn", "140", "A1133", 1965, "Physical Review", None, "Self-Consistent Equations Including Exchange and Correlation Effects"),
    (28, "Kotochigova", "55", "191", 1997, "Physical Review A", None, "Local-density-functional calculations of the energy of atoms"),
    (29, "Chakravorty", "47", "3649", 1993, "Physical Review A", None, "Ground-state correlation energies for atomic ions with 3 to 18 electrons"),
    (30, "Rodrigues", "86", "117", 2004, "Atomic Data and Nuclear Data Tables", "10.1016/j.adt.2003.11.005", None),
    (31, "Yerokhin", "44", "033103", 2015, "Journal of Physical and Chemical Reference Data", None, "Lamb shift of n=1 and n=2 states of hydrogen-like atoms, 1 <= Z <= 110"),
]
UA = {"User-Agent": "ref-check/1.0 (mailto:mawais9171@gmail.com)"}


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=40) as r:
        return json.load(r)


def fetch(doi, title, year):
    if doi:
        return get("https://api.crossref.org/works/" + urllib.parse.quote(doi))["message"]
    q = urllib.parse.quote(title)
    items = get(f"https://api.crossref.org/works?rows=5&query.bibliographic={q}"
                f"&filter=from-pub-date:{year},until-pub-date:{year}")["message"]["items"]
    return items[0] if items else None


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower().replace("ü", "u").replace("é", "e").replace("ó", "o"))


def main():
    bad = 0
    for n, fam, vol, page, year, jour, doi, title in REFS:
        try:
            it = fetch(doi, title, year)
        except Exception as e:
            print(f"[{n}] ERROR {e}")
            bad += 1
            continue
        if not it:
            print(f"[{n}] NOT FOUND in Crossref")
            bad += 1
            continue
        au = [a.get("family", "") for a in it.get("author", [])]
        yr = (it.get("issued", {}).get("date-parts", [[None]])[0] or [None])[0]
        pg = (it.get("page") or "").split("-")[0]
        cv = it.get("volume")
        cj = (it.get("container-title") or [""])[0]
        flags = []
        if au and norm(fam) not in norm(au[0]):
            flags.append(f"author:{au[0]}")
        if vol and cv and vol != cv:
            flags.append(f"volume:{cv}")
        if page and pg and page.lstrip("A") != pg.lstrip("A"):
            flags.append(f"page:{pg}")
        if yr != year:
            flags.append(f"year:{yr}")
        if jour and norm(jour.split()[0]) not in norm(cj + " ".join(it.get("container-title") or [])):
            flags.append(f"journal:{cj}")
        t = (it.get("title") or [""])[0]
        if title and norm(title)[:30] != norm(t)[:30]:
            flags.append("title-differs")
        status = "OK " if not flags else "CHECK"
        bad += bool(flags)
        print(f"[{n:2d}] {status} {fam} v{vol} p{page} {year} | Crossref: {', '.join(au[:3])}{'…' if len(au) > 3 else ''} | "
              f"{re.sub('<[^>]+>', '', t)[:70]} | {cj} {cv} {it.get('page')} {yr} | doi {it.get('DOI')}"
              + (f"  FLAGS {flags}" if flags else ""))
    print("\nreferences needing a look:", bad)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
