### Overview
This resource provides genome-wide maps of transcription initiation-derived cis-regulatory elements (CREs) identified using nucCAGE and the PRIME computational framework. The dataset includes cell line–specific predictions, FANTOM5 facet-level predictions, merged facet annotations, pooled CRE datasets, and TPM-normalized CAGE signal tracks.

CRE annotations are provided as coordinate-sorted BED9 files (bgzip-compressed and tabix-indexed) for efficient genomic querying and color-coded visualization in IGV and UCSC Genome Browser. Strand-specific bigWig files provide genome browser–ready visualization of transcription initiation activity across FANTOM5 facets and selected cell lines.

---

### Contents

#### PRIME_cellLines
CRE predictions for individual cell lines.

- Files are generated from genome-wide prediction outputs and filtered at PRIME score ≥ 0.75
- Each file corresponds to one cell line (e.g. GM12878, K562) and sample preparation approach, nuclei or whole cell, (e.g., K562_N, K562_C)
- Files are in BED9 format with PRIME scores scaled to 0–1000 and per-feature RGB coloring
- Each file contains a UCSC track line (`track ... itemRgb="On"`) for direct loading into IGV or UCSC Genome Browser with color-by-score support

---

#### PRIME_FANTOM5_facets
CRE predictions for individual FANTOM5 facets.

- Each file corresponds to a specific FANTOM5 cell, tissue, or cell line facet
- Files are in BED9 format with activity scores scaled to 0–1000 and per-feature RGB coloring
- Each file contains a UCSC track line for color-coded visualization
- File names follow the pattern:
  PRIME_FANTOM5_<facet>_0.5.bed.gz

---

#### PRIME_FANTOM5_agnostic
Pooled and merged CRE datasets across FANTOM5 facets.

Includes:
- Pooled CRE predictions across all facets
- Facet-merged CREs (aggregating activity across facets)
- CRE sets filtered at PRIME thresholds (0.5 and 0.75)
- Separate files for proximal, distal, and combined CREs

Merged facet files:
- Are in BED9+extra format: standard BED9 columns followed by a facets column
- Retain the maximum score across contributing facets
- Include facet annotations as a semicolon-separated list in the final column
- Contain a UCSC track line for color-coded visualization

Also includes:
- facet_statistics.tsv: summary statistics across facets

---

#### PRIME_FANTOM5_TPM_bw
TPM-normalized CAGE tracks for all FANTOM5 facets.

- Each facet is represented by two strand-specific bigWig files (.plus.bw and .minus.bw)
- Signals are derived from CTSS-level quantification and normalized to TPM
- These files are intended for direct visualization in genome browsers (e.g. IGV, UCSC Genome Browser)

---

#### PRIME_cellLines_TPM_bw
Pooled TPM-normalized CAGE signal tracks for individual cell lines.

- Cell lines: GM12878, K562, HCT116, HepG2, A549
- Each pooled group is represented by two strand-specific bigWig files (.plus.bw and .minus.bw)
- Pooling is defined by the **Type** column in each design matrix under `1.CTSSs/`
- For cell lines with multiple input amounts (e.g. GM12878 nuclei at 500K, 1M, 5M, 10M cells), each amount is kept as a separate pooled track
- Produced by `14.Resource/cellline_bigwigs_TPM.R` using saved CTSS RDS objects from `1.CTSSs/1.1.CTSSs.processing.Rmd` and `PRIME::writeBw`

---

### File format

All BED files follow standard BED conventions.

#### General properties
- Genome build: GRCh38
- Chromosome naming: UCSC style (chr1–chr22, chrX)
- Coordinate system: 0-based, half-open
- Sorting: sorted by chromosome and start coordinate
- Compression: bgzip
- Indexing: tabix (preset BED; use -p bed)

---

#### BED9 format (cell line and facet files)

Columns:
1. chrom
2. start
3. end
4. PRIME_id (chr:start-end)
5. score (integer, 0–1000)
6. strand (*)
7. thickStart (= start)
8. thickEnd (= end)
9. itemRgb (R,G,B — red scale; see Score definition)

Each file begins with a UCSC track line:

    track type=bed name="..." description="PRIME CREs (...)" visibility=2 itemRgb="On"

This line is skipped by tabix and enables per-feature coloring when the file is loaded as a custom track in IGV or UCSC Genome Browser.

---

#### Extended BED9 format (merged facet files)

Columns 1–9 as above, plus:
10. facets (semicolon-separated FANTOM5 facet labels)

---

### Score definition

The score column represents PRIME score scaled to the range 0–1000:

- Raw PRIME score values are bounded to [0,1]
- Scores are computed as:
  score = int(activity × 1000)

Interpretation:
- 0 = no confidence
- 1000 = maximal confidence

Thresholds used in the manuscript:
- 0.5 → score ≥ 500
- 0.75 → score ≥ 750

---

### Color encoding (itemRgb)

Each CRE feature is colored using a **red scale** that maps the PRIME score to an RGB value:

- score   0 → RGB `255,255,255` (white)
- score 1000 → RGB `255,0,0`   (red)
- Linear mapping: `R = 255`, `G = B = round(255 × (1 − score/1000))`

Higher-confidence CREs appear more red; lower-confidence features appear lighter/white. This encoding is active when files are loaded with the embedded UCSC track line (`itemRgb="On"`).

---

### Usage

Example query using tabix:

    tabix PRIME_cellLines/PRIME_GM12878_N_5M_0.75.bed.gz chr1:1000000-1100000

Files can be directly visualised in genome browsers such as IGV or UCSC Genome Browser. When loaded as a custom track in UCSC or via the track line in IGV, features are automatically colored by PRIME score (red scale).

---

### License

Data are released under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

---

### Citation

If you use this resource, please cite:
1. Einarsson, Navamajiti, et al, Mapping active regulatory elements from transcription initiation events, 2026, BioRxiv
2. This Zenodo record (doi: 10.5281/zenodo.19331480)