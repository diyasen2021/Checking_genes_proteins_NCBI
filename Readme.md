# Gene Table Verification — Setup & Run

## 1. Install dependencies
```bash
pip install biopython requests openpyxl pandas
```

## 2. Edit the config block at the top of `verify_gene_table.py`
- `INPUT_XLSX` — path to your `Genes_10th.xlsx`
- `NCBI_EMAIL` — **required** by NCBI's usage policy, put your real email
- `NCBI_API_KEY` — optional, free from your NCBI account, raises the rate
  limit from 3→10 requests/sec (58 rows × 2 NCBI calls each will take
  ~40s without a key, ~13s with one)

## 3. Run it
```bash
python verify_gene_table.py
```

Takes a few minutes (rate-limited to be polite to NCBI/UniProt). Produces
`gene_table_verification_report.xlsx` with three tabs:
- **Table1_Original** — your data, untouched
- **All_Checks** — every check run, per record, with PASS/FAIL/WARN/SKIP
- **Issues_Only** — just the rows worth looking at

## What it checks
**Offline** (no network): CDS/genomic sequence lengths vs stated lengths,
CDS→protein translation correctness, protein length, non-standard
characters, molecular weight recomputed from sequence, domain
start/end sanity, and stray whitespace in ID fields.

**Online** (needs internet): NCBI Nucleotide lookup for `NCBI_Accession`
+ sequence match, NCBI Protein lookup for `Protein_accession` + sequence
match, a best-effort InterPro domain check via UniProt's RefSeq
cross-reference, and a check that the NCBI record's title mentions the
stated `Gene_Name`/`Abbreviation` and that its organism matches `Species`
(catches an accession that resolves and even has a matching sequence,
but actually belongs to the wrong gene).

**Not automatically checked:** `Gene_Locus_ID`. Your table mixes AT-locus
tags, LOC ids, cultivar-specific codes, and free-text identifiers across
4 species — there's no single database endpoint that reliably resolves
all of those formats, so this field isn't spot-checked automatically.
The NCBI accession checks are the reliable substitute.

## Findings already surfaced (from the offline pre-run)
Ran the no-network checks against your file already — 18 flags out of 58
records, worth a look before you even run the online pass:

- **Likely real data errors (MW mismatches >3%):**
  - Record 4, **WAK2**: stated 42.32 kDa vs. recalculated 81.64 kDa for a
    732-aa protein — the recalculated value fits the protein length; the
    stated MW looks like a copy-paste error from a different record.
  - Record 3, **WAK1**: stated 95.83 kDa vs. recalculated 81.21 kDa (15% off)
  - Record 1, **GRP-3**: stated 14.89 vs. 14.29 kDa (4% off)
  - Record 42, **ERF1**: stated 23.87 vs. 24.69 kDa (3.5% off)
- **CDS doesn't translate cleanly to the stated protein:**
  - Record 14, **GLR2.8**: translated CDS is 948 aa, stated protein is 947 aa,
    first mismatch at position 1 — worth checking for an off-by-one or a
    stray residue.
- **CDS not found as a substring of the genomic/full sequence** (may be fine
  if the full sequence given is a different isoform/region, but worth a
  glance): Record 6 (**Stb6**), Record 33 (**WRKY33**)
- **Stray whitespace/newlines in ID fields** (cosmetic, but will break exact
  lookups if not cleaned first): records 5, 8, 11, 15, 17, 27, 36, 40, 55, 59
  — full detail in the CSV.

Full detail (field, issue type, description) is in
`offline_check_report_prerun.csv`.
