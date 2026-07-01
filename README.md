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

This is a lightweight **glue repository**: it holds the shared data and analysis directly and links to the annotation engines, but it does not contain their code. Each engine lives in its own independent GitHub repository, listed in [`engines.json`](engines.json) and pulled in by [`scripts/clone-all.sh`](scripts/clone-all.sh) for local development. They share data through the filesystem (genomes and the species table under `data/`) rather than by importing each other's code.

Engine repositories:

- [isoquant_annotator](https://github.com/Janek21/isoquant_annotator) — IsoQuant long-read assembly + annotation guided
- [LyRic_annotator](https://github.com/Janek21/LyRic_annotator) — LyRic wrapper: long-read assembly + annotation unguided
- [geneid-training](https://github.com/Janek21/geneid-training) — geneid train, predict, and AGAT merge

```
parent/
  README.md                this file
  engines.json               manifest of the engine repositories
  scripts/clone-all.sh     clones/updates every engine into this directory
  scripts/tables_setup.sh   obtains all current available protist data by strain and saves in data/
  data/                   acquire long-read RNA-seq (ENA) + genomes (NCBI); species tables
  busco_references/        BUSCO baselines for reference genomes + shared taxonomy utilities
  protist_table/           protist species survey
  result_analytics/        aggregate per-engine metrics, compare, and report
  isoquant_annotator/      cloned engine (not tracked here)
  LyRic_annotator/         cloned engine (not tracked here)
  geneid-training/         cloned engine (not tracked here)
```

The shared directories (`data/`, `busco_references/`, `protist_table/`, `result_analytics/`) are part of this repository. The engine directories are cloned on demand and ignored by Git, so they stay independent and the glue repo only records which commit each engine should be at, via the manifest.

## Setup

### Engine setup

```bash
git clone https://github.com/Janek21/<this-repo>.git protist_annotation
cd protist_annotation
bash scripts/clone-all.sh
```

`clone-all.sh` reads `engines.json` and clones each engine into `protist_annotation/` at the top level, or updates it if already present.

### Species setup

```bash
bash scripts/tables_setup.sh
cd data

#download species
sbatch scripts/referenceDownload.sh

#rename directory taxids to match the table
python scripts/rename_to_longread_taxid.py --species-dir species --longread longread_protists.tsv --apply
```

All strains that currently present an assembly and long-read RNA-seq data are downloaded to `data/species`.

#### Evaluation of gathered species

Gather the list of downloaded species names(without the taxid):
```bash
ls data/species/| sed 's/_[0-9]\+$//' > busco_references/dataspecie.txt
cd busco_references
```

Adapt `busco_references/turbo.sh` to the list of species; the array number has to reflect the amount of species in busco_references/dataspecie.txt
```bash
wc -l dataspecie.txt

#modify turbo.sh array number
sbatch turbo.sh
```
_If the goal is to evaluate only 1 specie use `busco_references/monoSpecie.sh` instead of `busco_references/turbo.sh`_


## Dataset

All available(65 currently) protist species with public long-read RNA-seq in the NCBI SRA (currently 354 accessions, Oxford Nanopore and PacBio), spanning diverse protist lineages. The current state of protist strain annotation can be visualized in the [strains](protist_table/strain_protists.html) and [accessions](protist_table/protists.html) tables.

## Annotation production

The annotation engines are designed to be run LyRic>Isoquant>geneid as Isoquant can reuse the raw data from LyRic and geneid is a final refination process.

For running each engine you can find the instructions in their respective README, some cases might require a species list, for those cases use the generated `busco_references/dataspecie.txt`.
```bash
cp busco_references/dataspecie.txt *annotator
```

## Result collection

Once all annotation engines have been run the results can be collected into one [general table](result_analytics/overview.html) by:

 ```bash
 python result_analytics/build_results_table.py --html --out result_analytics/overview.tsv

 #the busco stats(.png) can be gathered with
 python result_analytics/collect_busco_summaries.py --refresh
 ```

## License

Code is under its respective repository licenses; the written report is under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 license. Copyright (c) Jan Izquierdo i Ramos.
