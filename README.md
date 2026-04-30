Code for analyses and figures in the paper "Mapping active cis-regulatory elements from transcription initiation events" by Einarsson, Navamajiti, et al. 2026

PRIME predictions in K562, GM12878, HepG2, HCT116, A549 cells and FANTOM5 facets as well as CAGE signal tracks are available at Zenodo (doi: [10.5281/zenodo.19712783](https://zenodo.org/records/19712783)).

The code in this repository relies on the **PRIME toolkit**, three interconnected tools for the analysis of
transcription initiation data (e.g., CAGE):

| Tool | Type | Purpose |
|---|---|---|
| [PRIMEprep](https://github.com/anderssonlab/PRIMEprep) | Bash pipeline | Raw FASTQ → QC → trimming → mapping → BigWig |
| [PRIME](https://github.com/anderssonlab/PRIME) | R package | CTSS quantification, divergent loci, promoter decomposition, normalization, noise estimation |
| [PRIMEmodel](https://github.com/anderssonlab/PRIMEmodel) | R package + Python | Genome-wide prediction of regulatory elements |
