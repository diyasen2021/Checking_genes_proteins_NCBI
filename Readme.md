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

