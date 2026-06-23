# Protist annotation: a methods improvement

A methodological improvement for protist gene annotation built around the **geneid** gene predictor, driven by long-read RNA-seq evidence. Master's thesis (TFM) in Data Science, developed at the Center for Genomic Regulation (CRG).

## Overview

Most ab initio predictors are trained on a few model organisms, so they annotate divergent protist lineages poorly (wrong codon usage, intron structure, and in some lineages non-standard genetic codes). This project combines long-read transcript evidence with species-specific geneid training in one workflow:

1. Assemble long-read RNA-seq into transcript models (LyRic and IsoQuant).
2. Select the best-supported annotation using BUSCO plus gene/transcript counts
   and genome size.
3. Use it as the base annotation (against the reference where one exists).
4. Train geneid on that base and run it across the genome.
5. Merge the geneid predictions back with the base using AGAT.

## How this repository is organized

This is a lightweight **glue repository**: it holds the shared data and analysis directly and links to the annotation engines, but it does not contain their code. Each engine lives in its own independent GitHub repository, listed in [`repos.json`](repos.json) and pulled in by [`scripts/clone-all.sh`](scripts/clone-all.sh) for local development. They share data through the filesystem (genomes and the species table under `data/`) rather than by importing each other's code.

Engine repositories:

- [isoquant_annotator](https://github.com/Janek21/isoquant_annotator) — IsoQuant long-read assembly + annotation
- [LyRic_annotator](https://github.com/Janek21/LyRic_annotator) — LyRic wrapper; clones the LyRic engine ([LyRic_nonhuman](https://github.com/Janek21/LyRic_nonhuman)) per run
- [geneid-training](https://github.com/Janek21/geneid-training) — geneid train, predict, and AGAT merge

```
repo/
  README.md              this file
  repos.json             manifest of the engine repositories
  scripts/clone-all.sh   clones/updates every engine into this directory
  data/                  acquire long-read RNA-seq (ENA) + genomes (NCBI); species tables
  busco_references/      BUSCO baselines for reference genomes + shared taxonomy utilities
  protist_table/         protist species survey
  result_analytics/      aggregate per-engine metrics, compare, and report
  isoquant_annotator/    cloned engine (not tracked here)
  LyRic_annotator/       cloned engine (not tracked here)
  geneid-training/       cloned engine (not tracked here)
```

The shared directories (`data/`, `busco_references/`, `protist_table/`, `result_analytics/`) are part of this repository. The engine directories are cloned on demand and ignored by Git, so they stay independent and the glue repo only records which commit each engine should be at, via the manifest.

## Setup

```bash
git clone https://github.com/Janek21/<this-repo>.git repo
cd repo
bash scripts/clone-all.sh
```

`clone-all.sh` reads `repos.json` and clones each engine into `repo/` at the top level, or updates it if already present.

## Dataset

64 protist species with public long-read RNA-seq in the NCBI SRA (354 accessions, Oxford Nanopore and PacBio), spanning diverse protist lineages.

## License

Code is under its respective repository licenses; the written report is under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 license. Copyright (c) Jan Izquierdo i Ramos.
