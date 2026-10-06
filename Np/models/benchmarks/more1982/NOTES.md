# More (1982) screened hydrogenic model: benchmark status

**Status: constants NOT OBTAINABLE within the time box (about 30 min, 2026-10-05).** No model was implemented
and nothing was scored. No constants were reconstructed from memory.

## Citation (verified)
R. M. More, "Electronic energy-levels in dense plasmas", J. Quant. Spectrosc. Radiat. Transfer **27**(3), 345-357
(March 1982), doi:10.1016/0022-4073(82)90127-3. Verified with the Crossref API record for that DOI (title, author,
journal, volume 27, issue 3, pages 345-357, 1982-03).

## What was found (model form only, no constants)
IPPJ-AM-31 (Kagawa, Kato, Watanabe, "Atomic processes in hot dense plasmas", Nagoya IPP report,
http://dpc.nifs.ac.jp/IPPJ-AM/IPPJ-AM-31.pdf, pp. 9-10, Eqs. 15-18) summarises the SHM as follows.
The orbital energy is E_n = E0_n - Q_n^2 e^2/(2 a0 n^2), with Q_n = Z - sum_{m<n} sigma(n,m) P_m - (1/2) sigma(n,n) P_n.
The outer-screening shift is E0_n = -(1/2) sigma(n,n) P_n e^2/r_n ... + sum_{m>n} P_m sigma(m,n) e^2/r_m,
with r_n = a0 n^2/Q_n. The source cites a 1981 LLNL report by More. It does **not** print the sigma(n,m) table.
The OCR of that page is noisy, so treat the equation transcription as indicative only.

## Search log (where we looked)
- ScienceDirect abstract page (0022407382901273): HTTP 403 (paywalled).
- OSTI.gov: every URL refused the connection from this machine (ECONNREFUSED / timeout). That covers the ETDEWEB
  conference record 6694951 "Electronic energy-levels in dense plasmas", UCRL-93926 "Atoms in dense plasmas"
  (purl 6278776), UCRL-88511 and ETDEWEB 10110158 "Screening constants for plasma". These are the most likely
  open copies of the table and could not be reached.
- NIFS repository record 9860 (NIFS-235, Kawata/Kato/Kiyokawa 1993, "Screening constants for plasma", which compares
  against More's constants): metadata only, no full text.
- IPPJ-AM-31 (NIFS): equations only, no table (see above).
- Semantic Scholar citations API for the DOI: 187 citing papers. We scanned the open-access or arXiv ones that might
  quote the table, all read in memory:
  - arXiv 0709.4473, 2007.05243, 2210.04938, 2211.16464 (SpK), 2401.03180 (MolDStruct), 2406.06233,
    1404.4531, 2103.07663, 2602.22064: cite More 1982, no table.
  - Mendoza et al. 2011 (oa.upm.es/11165, INVE_MEM_2011_102088.pdf): compares with More's SHM but prints only its
    own constants.
  - Di Rocco/Lanzini, Braz. J. Phys. 46 (2016) 175 (Redalyc): no More table.
  - SciELO Braz. J. Phys. 2005, "An efficient screening approach…": cites More, no table.
- French thesis 2017SACLS543 (theses.hal.science/tel-01778583): blocked by an anti-bot JavaScript challenge, which we
  did not bypass.
- ULPGC accedacris records 10553/45566, 45552, 46068, 45575 (Rubiano/Mendoza group): HTTP 403 or timeout.
- CONICET repository 11336/115941: connection refused.
- UNT digital library (mirrors many LLNL reports): the search endpoint returned nothing.
- Generic web searches (More screening constants table, sigma(n,m), GitHub implementations): no verbatim, verifiable
  reproduction of the table.

## How to finish this later
Get the table from the original paper (library access) or from an open LLNL preprint on OSTI, using a network that
can reach osti.gov. Transcribe it into `constants.json` with the page and table number. Then validate against the
example numbers printed in the paper before scoring.
