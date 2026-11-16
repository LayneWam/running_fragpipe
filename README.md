# AAS Stabilization: FragPipe Reprocessing Pipeline

Pipeline for downloading TMT proteomics raw data from the [Proteomic Data Commons (PDC)](https://pdc.cancer.gov/), organizing it by plex, and running per-plex [FragPipe](https://github.com/Nesvilab/FragPipe) searches on a Slurm cluster.

## Table of Contents

- [Directory Layout](#directory-layout)
  - [Scratch](#scratch-scratchlamwayn)
  - [Home](#home-homelamwaynaas_stabilizationreprocessing)
- [Software Setup](#software-setup)
- [Workflow](#workflow)
- [Usage](#usage)
  - [Default (per-plex FASTA)](#default-per-plex-fasta)
  - [Single-FASTA mode](#single-fasta-mode)
  - [Testing a single plex](#testing-a-single-plex)
- [Notes](#notes)

---

## Directory Layout

### Scratch: `/scratch/lam.wayn`

Contains raw data, results, FASTA files, metadata, and logs.

```
/scratch/lam.wayn/
├── FASTA/                                    # FASTA files
│   ├── per_plex/                             # Per-plex FASTA files
│   ├── per_plex_fragpipe/                    # FASTA files with decoys
│   ├── list_mutated_seqs_beneficialMuts.fasta  # Beneficial mutations
│   └── multi_singleton_ref_supp7.fasta       # Singles, multis, and refs
├── logs/                                     # Job logs
├── METADATA/                                 # Downloaded manually from PDC Commons
│   └── <experiment>/
│       ├── PDC file manifest .csv            # Raw file → plex mapping
│       ├── PDC study biospecimen .csv        # Sample list
│       └── PDC study experimental .csv       # Plex → channel → sample mapping
├── MS_SEARCH/
│   ├── annotation/                           # Plex annotation .tsv (sample, plex, channel, reference)
│   ├── manifests/                            # Raw file → plex mapping
│   ├── results/
│   │   └── <plex>/                           # Per-plex results
│   ├── workflow/                             # Per-plex workflow files
│   ├── fragpipe.workflow                     # Sample workflow file
│   ├── plex_list.txt                         # List of plexes
│   └── search_setup_summary                  # Summary of raw file search (how many were found)
└── RAW/
    └── <plex>/
        ├── *annotation.txt                   # Channel → sample mapping and reference identification
        └── *.mzML, *.mzBIN                   # Raw data
```

### Home: `/home/lam.wayn/aas_stabilization/reprocessing`

Contains software, raw data download tools, and scripts.

```
/home/lam.wayn/aas_stabilization/reprocessing/
├── bin/
│   ├── fragpipe/                     # FragPipe (Linux) + MSFragger, IonQuant, diaTracer
│   └── jdk/                          # Java (Temurin JDK 17, x64 Linux HotSpot)
├── download/
│   ├── download_files.sh             # Downloads raw files using the consolidated file manifest
│   └── PDC_meta/                     # Consolidated metadata
│       ├── pdc_all_files.tsv         # Manifest of all raw files across all experiments
│       ├── pdc_file_tmt_map.tsv      # Raw file → plex and channel → sample, across all experiments
│       └── pdc_meta.tsv              # Sample manifest across all experiments
└── SCRIPTS/
    ├── consolidate_metadata.py       # Consolidates metadata from all experiments (used to download raw files)
    ├── group_raw.py                  # Groups raw files into plex subdirectories
    ├── run_ms_search.sh              # Master/main script
    ├── add_decoys.py                 # Adds decoys to FASTA, saves to FASTA/per_plex_fragpipe/
    ├── generate_manifests.py         # Generates per-plex workflows, annotations, manifests, and plex_list
    └── submit_fragpipe.sh            # Runs FragPipe (set JAVA_HOME, FRAGPIPE_BIN, FRAGPIPE_TOOLS here)
```

---

## Software Setup

Install the following into `bin/`:

| Software | Version | Source |
|---|---|---|
| FragPipe (Linux) | latest | [GitHub releases](https://github.com/Nesvilab/FragPipe/releases) |
| MSFragger | 4.4.1 | [msfragger.arsci.com/msfragger](https://msfragger.arsci.com/msfragger/) |
| IonQuant | 1.11.20 | [msfragger.arsci.com/ionquant](https://msfragger.arsci.com/ionquant/) |
| diaTracer | 2.2.1 | [msfragger-upgrader.nesvilab.org/diatracer](https://msfragger-upgrader.nesvilab.org/diatracer/) |
| Java JDK (x64_linux_hotspot) | 17.0.21+9 | [Adoptium Temurin 17](https://github.com/adoptium/temurin17-binaries/releases/tag/jdk-17.0.21+9) |

After installing, set the following variables in `SCRIPTS/submit_fragpipe.sh`:

```bash
JAVA_HOME=...
FRAGPIPE_BIN=...
FRAGPIPE_TOOLS=...
```

---

## Workflow

1. **Download software:** FragPipe (with MSFragger, IonQuant, diaTracer) and Java. See [Software Setup](#software-setup).
2. **Download metadata manually** from PDC Commons into `/scratch/lam.wayn/METADATA/<experiment>/`.
3. **Consolidate metadata** to create the combined file manifest:
   ```bash
   python consolidate_metadata.py
   ```
4. **Download raw data** programmatically:
   ```bash
   bash download_files.sh
   ```
5. **Group raw data** into plex subdirectories:
   ```bash
   python group_raw.py
   ```
6. **Run the MS search:**
   ```bash
   bash run_ms_search.sh
   ```
   This runs, in order:
   1. `add_decoys.py`
   2. `generate_manifests.py`
   3. `submit_fragpipe.sh`

---

## Usage

All commands below are run from `/home/lam.wayn/aas_stabilization/reprocessing/SCRIPTS`.

### Default (per-plex FASTA)

By default, each plex is searched against its own FASTA, since each plex contains different samples with different protein variants.

```bash
bash run_ms_search.sh
```

Jobs are submitted as a Slurm array, 5 at a time (Slurm limit):

```bash
sbatch --array=1-{NUM_PLEXES}%5 submit_fragpipe.sh
```

### Single-FASTA mode

To compare all plexes against a single common FASTA:

```bash
bash run_ms_search.sh --single-fasta /path/to/your.fasta /path/to/your.workflow
```

Example:

```bash
bash run_ms_search.sh --single-fasta \
    /scratch/lam.wayn/FASTA/list_mutated_seqs_beneficialMuts.fasta \
    /scratch/lam.wayn/MS_SEARCH/fragpipe.workflow
```

### Testing a single plex

To test-run only the final FragPipe step on one plex, pass the index of that plex (from `plex_list.txt`):

```bash
sbatch --array={PLEX_INDEX} submit_fragpipe.sh
```

---

## Notes

- **PDC manifest files expire after 7 days.** Re-download them after that.
- **Annotation files** must be placed in the plex raw data folder alongside the raw data so FragPipe knows which channel is the reference, and their names must end with `annotation.txt`.
- **Each annotation file can only have a single reference**, set by `ref_tag` in the workflow file.
- **TMT type** is set automatically in the workflow file by `generate_manifests.py`.
- Set **`add_Ref=-1`** in the workflow to use the real reference channel instead of a virtual reference.
