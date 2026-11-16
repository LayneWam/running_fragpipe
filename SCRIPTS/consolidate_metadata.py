#!/usr/bin/env python3
"""
Consolidate PDC metadata from all tissue folders.

Each tissue folder (TISSUE_DIR/<tissue>/) may contain one or more of:
- PDC_file_manifest_*.csv       : raw file download info (16 columns, see FILE_MANIFEST_COLUMNS)
- PDC_study_biospecimen_*.csv   : sample/patient metadata
- PDC_study_experimental_*.csv  : TMT/iTRAQ plex layout (one row per Run Metadata ID)

If a folder holds several exports of the same kind (e.g. an old and a new
download, or several studies), all are read and rows are de-duplicated,
keeping the row from the NEWEST file. This matters for the file manifest,
because only the newest export has unexpired download links.

Outputs (in META_DIR):
- pdc_all_files.tsv        : the 16 original PDC file-manifest columns, in PDC's
                             order, plus tissue_folder as column 17.
                             File Download Link is column 16 (1-based).
- pdc_meta.tsv             : one row per aliquot
- pdc_file_tmt_map.tsv     : raw file -> TMT channel -> sample

Usage:
    python consolidate_metadata.py
"""

import base64
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import parse_qs, urlparse

META_DIR = Path('/home/lam.wayn/aas_stabilization/reprocessing/PDC_meta')
TISSUE_DIR = Path('/scratch/lam.wayn/METADATA')

OUTPUT_FILE_MANIFEST = META_DIR / 'pdc_all_files.tsv'
OUTPUT_SAMPLE_META = META_DIR / 'pdc_meta.tsv'
OUTPUT_FILE_TMT_MAP = META_DIR / 'pdc_file_tmt_map.tsv'

# Exact header of PDC_file_manifest_*.csv exports (0-based index shown)
FILE_MANIFEST_COLUMNS = [
    'File ID',               # 0
    'File Name',             # 1
    'Run Metadata ID',       # 2
    'Protocol',              # 3
    'Study Name',            # 4
    'PDC Study ID',          # 5
    'PDC Study Version',     # 6
    'Study ID',              # 7
    'Project Name',          # 8
    'Data Category',         # 9
    'File Type',             # 10
    'Access',                # 11
    'File Size (in bytes)',  # 12
    'Md5sum',                # 13
    'Downloadable',          # 14
    'File Download Link',    # 15
]

# Columns we rely on in the other two exports
BIOSPECIMEN_REQUIRED = ['Aliquot ID', 'Aliquot Submitter ID', 'Sample Submitter ID',
                        'Case Submitter ID', 'Sample Type', 'Primary Site',
                        'Disease Type', 'Tissue Type']
EXPERIMENTAL_REQUIRED = ['Run Metadata ID']

# Plex channel columns in the experimental export start with one of these.
# Channel names vary by plex (TMT-10 uses 'tmt_131', TMT-11 adds 'tmt_131c',
# TMT-16/18 add 'tmt_132n'... ), so they are detected from the header.
CHANNEL_PREFIXES = ('tmt_', 'itraq_')


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------

def find_csv_files(folder, prefix):
    """All CSVs in folder starting with prefix, newest first."""
    return sorted(folder.glob(f'{prefix}*.csv'),
                  key=lambda p: p.stat().st_mtime, reverse=True)


def read_csv_with_bom(filepath):
    """Read a PDC CSV export (UTF-8 with BOM). Returns (header, rows)."""
    with open(filepath, 'r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return reader.fieldnames or [], rows


def check_columns(filepath, header, required):
    missing = [c for c in required if c not in header]
    if missing:
        print(f"  WARNING: {filepath.name} is missing columns: {missing}")
        return False
    return True


def read_exports(tissue_dir, prefix, required, key_col):
    """Read every export of one kind in a folder, newest first, de-duplicated on key_col."""
    rows_out, seen = [], set()
    files = find_csv_files(tissue_dir, prefix)
    if len(files) > 1:
        print(f"  NOTE: {len(files)} {prefix} files; newest wins on duplicates "
              f"({files[0].name})")
    for path in files:
        header, rows = read_csv_with_bom(path)
        if not check_columns(path, header, required):
            continue
        for r in rows:
            key = r.get(key_col, '')
            if key in seen:
                continue
            seen.add(key)
            r['tissue_folder'] = tissue_dir.name
            r['_source'] = path.name
            r['_header'] = header
            rows_out.append(r)
    return rows_out


def consolidate_all_tissues():
    all_files, all_biospecimens, all_experiments = [], [], []

    if not TISSUE_DIR.is_dir():
        sys.exit(f"Tissue directory not found: {TISSUE_DIR}")

    tissue_folders = sorted(d for d in TISSUE_DIR.iterdir() if d.is_dir())

    for tissue_dir in tissue_folders:
        print(f"Processing {tissue_dir.name}...")

        files = read_exports(tissue_dir, 'PDC_file_manifest',
                             FILE_MANIFEST_COLUMNS, 'File ID')
        n_before = len(files)
        files = [f for f in files if f.get('Downloadable', '').strip().lower() == 'yes'
                 and f.get('File Download Link', '').strip()]
        if n_before != len(files):
            print(f"  - skipped {n_before - len(files)} non-downloadable files")
        all_files.extend(files)
        print(f"  - {len(files)} files")

        samples = read_exports(tissue_dir, 'PDC_study_biospecimen',
                               BIOSPECIMEN_REQUIRED, 'Aliquot ID')
        all_biospecimens.extend(samples)
        print(f"  - {len(samples)} biospecimens")

        exps = read_exports(tissue_dir, 'PDC_study_experimental',
                            EXPERIMENTAL_REQUIRED, 'Run Metadata ID')
        all_experiments.extend(exps)
        print(f"  - {len(exps)} experimental runs")

    # A file can appear under more than one tissue folder; keep the first.
    deduped, seen = [], set()
    for f in all_files:
        if f['File ID'] in seen:
            continue
        seen.add(f['File ID'])
        deduped.append(f)
    if len(deduped) != len(all_files):
        print(f"\nRemoved {len(all_files) - len(deduped)} files listed in more than one tissue folder")

    return deduped, all_biospecimens, all_experiments


# ----------------------------------------------------------------------------
# Mapping
# ----------------------------------------------------------------------------

def build_case_sample_map(biospecimens):
    """Map (case_submitter_id, sample_type) -> list of aliquot records."""
    case_map = {}
    for b in biospecimens:
        case_id = b.get('Case Submitter ID', '').strip()
        if not case_id:
            continue
        key = (case_id, b.get('Sample Type', '').strip())
        case_map.setdefault(key, []).append({
            'aliquot_submitter_id': b.get('Aliquot Submitter ID', ''),
            'sample_submitter_id': b.get('Sample Submitter ID', ''),
            'tissue_type': b.get('Tissue Type', ''),
            'primary_site': b.get('Primary Site', ''),
        })
    return case_map


def parse_channel_cell(value):
    """Cells look like 'CaseID\\r\\nSample Type'; return (case_id, sample_type)."""
    lines = [l.strip() for l in value.replace('\r', '').split('\n') if l.strip()]
    case_id = lines[0] if lines else ''
    sample_type = lines[1] if len(lines) > 1 else ''
    return case_id, sample_type


def build_run_to_samples_map(experiments, case_map):
    """Map Run Metadata ID -> list of channel/sample records."""
    run_map = {}
    ambiguous = 0
    for exp in experiments:
        run_id = exp.get('Run Metadata ID', '').strip()
        if not run_id:
            continue
        channels = [c for c in exp['_header'] if c.lower().startswith(CHANNEL_PREFIXES)]

        samples_in_run = []
        for channel in channels:
            value = exp.get(channel, '')
            if not value or not value.strip():
                continue
            case_id, sample_type = parse_channel_cell(value)
            matches = case_map.get((case_id, sample_type), [])
            if len(matches) > 1:
                ambiguous += 1
            info = matches[0] if len(matches) == 1 else {}
            samples_in_run.append({
                'case_submitter_id': case_id,
                'sample_type': sample_type,
                'tmt_channel': channel,
                # Left blank when a case has several aliquots of this sample type,
                # because the experimental export doesn't say which one was used.
                'aliquot_submitter_id': info.get('aliquot_submitter_id', ''),
                'tissue_type': info.get('tissue_type', ''),
                'primary_site': info.get('primary_site', ''),
            })
        run_map[run_id] = samples_in_run

    if ambiguous:
        print(f"NOTE: {ambiguous} channel entries match several aliquots; "
              f"aliquot_submitter_id left blank for those")
    return run_map


# ----------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------

def link_expiry(url):
    """Expiry time of a CloudFront signed URL, or None if it can't be read."""
    try:
        policy = parse_qs(urlparse(url).query)['Policy'][0]
        policy = policy.replace('-', '+').replace('_', '=').replace('~', '/')
        data = json.loads(base64.b64decode(policy + '==='))
        epoch = data['Statement'][0]['Condition']['DateLessThan']['AWS:EpochTime']
        return datetime.fromtimestamp(epoch, timezone.utc)
    except Exception:
        return None


def write_file_manifest(files):
    """Original 16 PDC columns in PDC order, then tissue_folder (column 17)."""
    headers = FILE_MANIFEST_COLUMNS + ['tissue_folder']

    with open(OUTPUT_FILE_MANIFEST, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(headers)
        for file in files:
            # Tabs/newlines inside a value would break the TSV
            writer.writerow([str(file.get(c, '')).replace('\t', ' ')
                             .replace('\r', ' ').replace('\n', ' ').strip()
                             for c in headers])

    print(f"\nWrote {len(files)} files to {OUTPUT_FILE_MANIFEST}")

    expiries = [e for e in (link_expiry(f['File Download Link']) for f in files) if e]
    if expiries:
        soonest = min(expiries)
        now = datetime.now(timezone.utc)
        status = "ALREADY EXPIRED" if soonest <= now else f"{soonest - now} from now"
        print(f"Earliest download-link expiry: {soonest:%Y-%m-%d %H:%M} UTC ({status})")
        if soonest <= now:
            print("  -> re-export file manifests from PDC before downloading.")


def write_sample_metadata(biospecimens):
    headers = ['aliquot_submitter_id', 'sample_submitter_id', 'case_submitter_id',
               'sample_type', 'tissue_type', 'primary_site', 'disease_type', 'tissue_folder']

    with open(OUTPUT_SAMPLE_META, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(headers)
        for b in biospecimens:
            writer.writerow([
                b.get('Aliquot Submitter ID', ''),
                b.get('Sample Submitter ID', ''),
                b.get('Case Submitter ID', ''),
                b.get('Sample Type', ''),
                b.get('Tissue Type', ''),
                b.get('Primary Site', ''),
                b.get('Disease Type', ''),
                b.get('tissue_folder', ''),
            ])

    print(f"Wrote {len(biospecimens)} samples to {OUTPUT_SAMPLE_META}")


def write_file_tmt_mapping(files, run_map):
    headers = ['file_name', 'run_metadata_id', 'tmt_channel', 'case_submitter_id',
               'sample_type', 'tissue_type', 'aliquot_submitter_id', 'tissue_folder']

    rows_written = 0
    unmapped = 0
    with open(OUTPUT_FILE_TMT_MAP, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(headers)
        for file in files:
            file_name = file.get('File Name', '')
            run_id = file.get('Run Metadata ID', '').strip()
            tissue_folder = file.get('tissue_folder', '')
            samples = run_map.get(run_id, [])

            if samples:
                for s in samples:
                    writer.writerow([file_name, run_id, s['tmt_channel'],
                                     s['case_submitter_id'], s['sample_type'],
                                     s['tissue_type'], s['aliquot_submitter_id'],
                                     tissue_folder])
                    rows_written += 1
            else:
                unmapped += 1
                writer.writerow([file_name, run_id, '', '', '', '', '', tissue_folder])
                rows_written += 1

    print(f"Wrote {rows_written} file-sample mappings to {OUTPUT_FILE_TMT_MAP}")
    if unmapped:
        print(f"NOTE: {unmapped} files have no matching run in any experimental export")


def main():
    print("Consolidating PDC metadata from all tissue folders...\n")

    all_files, all_biospecimens, all_experiments = consolidate_all_tissues()

    print("\nTotals:")
    print(f"  Files: {len(all_files)}")
    print(f"  Biospecimens: {len(all_biospecimens)}")
    print(f"  Experimental runs: {len(all_experiments)}")

    case_map = build_case_sample_map(all_biospecimens)
    run_map = build_run_to_samples_map(all_experiments, case_map)

    META_DIR.mkdir(parents=True, exist_ok=True)
    write_file_manifest(all_files)
    write_sample_metadata(all_biospecimens)
    write_file_tmt_mapping(all_files, run_map)

    print("\nDone!")


if __name__ == '__main__':
    main()
