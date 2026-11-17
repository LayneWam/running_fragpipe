# -*- coding: utf-8 -*-
"""
generate_manifests.py

Generate FragPipe manifests and TMT annotations for per-plex MS searches.

Each TMT plex is searched independently against its custom FASTA database
(reference proteome + plex-specific mutant + compensatory entries).

Reads the TMT mapping file to group RAW files by plex, then generates:
  - Per-plex FragPipe manifest (.fp-manifest)
  - Per-plex TMT channel annotation, saved in the plex's RAW subdirectory
    (RAW_DIR/<plex>/<plex>_annotation.txt)
  - Per-plex workflow (from template, with FASTA path patched)
  - Plex list for SLURM array submission

Usage:
    python3 generate_manifests.py \
        --tmt-map /path/to/pdc_file_tmt_map.tsv \
        --raw-dir /path/to/RAW/ \
        --fasta-dir /path/to/FASTA/per_plex/ \
        --workflow-template /path/to/template.workflow \
        -o /path/to/MS_SEARCH/

    # Single mode: search every plex against one shared FASTA instead
    python3 generate_manifests.py \
        --tmt-map /path/to/pdc_file_tmt_map.tsv \
        --raw-dir /path/to/RAW/ \
        --single-fasta /path/to/combined_decoy.fasta \
        --workflow-template /path/to/template.workflow \
        -o /path/to/MS_SEARCH/
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


# TMT channel name mapping: pdc_file_tmt_map format -> FragPipe format
TMT_CHANNEL_MAP = {
    "tmt_126": "126",
    "tmt_127n": "127N",
    "tmt_127c": "127C",
    "tmt_128n": "128N",
    "tmt_128c": "128C",
    "tmt_129n": "129N",
    "tmt_129c": "129C",
    "tmt_130n": "130N",
    "tmt_130c": "130C",
    "tmt_131" : "131N",
    "tmt_131n": "131N",
    "tmt_131c": "131C"
}

def detect_tmt_type(channel_map):
    return "TMT-11" if "tmt_131c" in channel_map else "TMT-10"

def load_tmt_map(tmt_map_path):
    """
    Load TMT mapping file.

    Returns:
        plex_to_files: dict of run_metadata_id -> set of file_names
        plex_to_channels: dict of run_metadata_id -> list of
            {"file_name", "channel", "case_id", "sample_type"}
    """
    plex_to_files = defaultdict(set)
    plex_to_channels = defaultdict(list)

    with open(tmt_map_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            run_id = row["run_metadata_id"].strip()
            file_name = row["file_name"].strip()
            channel = row.get("tmt_channel", "").strip().lower()
            case_id = row.get("case_submitter_id", "").strip()
            sample_type = row.get("sample_type", "").strip()

            plex_to_files[run_id].add(file_name)
            plex_to_channels[run_id].append({
                "file_name": file_name,
                "channel": channel,
                "case_id": case_id,
                "sample_type": sample_type,
            })

    return plex_to_files, plex_to_channels

def _is_reference_case(case_id):
    cid = case_id.lower()
    return cid == any(tag in cid for tag in ("pool", "ref", "global"))

def _normalize_channel(ch):
    c = str(ch).strip().lower()
    if c.startswith("tmt_"):
        c = c[len("tmt_"):]
    return c


def get_unique_channel_annotation(channel_entries, plex_id):
    """
    Collapse channel entries to unique (channel -> sample) mapping.

    Multiple file_names (fractions) map to the same TMT channels, so
    deduplicate to one row per channel for the annotation file.
    """
    channel_to_sample = {}
    ref_cand = []
    for entry in channel_entries:
        ch = entry["channel"]
        if ch not in channel_to_sample:
            case_id = entry["case_id"]
            sample_type = entry["sample_type"]
            cid = case_id.lower()
            if cid == "" or any(tag in cid for tag in ("ref", "pool", "global")):
                ref_cand.append(ch)
            label = f"{case_id}_{sample_type}" if sample_type else case_id
            condition = sample_type if sample_type else "unknown"
            channel_to_sample[ch] = {
                "sample": label,
                "condition": condition,
            }

    reference = {"sample": "Reference", "condition": "reference"}
    if len(ref_cand) < 1:
        print(f"No reference found for {plex_id}")
    elif len(ref_cand) == 1:
        channel_to_sample[ref_cand[0]] = reference
    else:
        ref = next((ch for p in ("131c", "126") for ch in ref_cand
                    if _normalize_channel(ch) == p), None)
        if ref is None:
            rank = {_normalize_channel(k): i for i, k in enumerate(TMT_CHANNEL_MAP)}
            ref = max(ref_cand, key=lambda ch: rank.get(_normalize_channel(ch), -1))
            print(f"annotation warning for {plex_id}")
        channel_to_sample[ref] = reference
    return channel_to_sample


def write_manifest(raw_files, plex_id, out_path):
    """Write FragPipe manifest file."""
    with open(out_path, "w") as f:
        for raw_path in sorted(raw_files):
            # experiment = plex_id, bioreplicate = 1, data_type = DDA
            f.write(f"{raw_path}\t{plex_id}\t1\tDDA\n")


def write_annotation(channel_map, plex_id, out_path, tmt_type="TMT-10"):
    """
    Write the TMT channel annotation file in the format FragPipe reads:
    two whitespace-separated columns, channel then sample name, no header.

        126   11BR047_Primary_Tumor
        127N  11BR043_Primary_Tumor
        ...
        131   Reference

    FragPipe splits each line on whitespace, so spaces inside sample names
    are replaced with underscores. TMT-10 has a single 131 channel, labelled
    "131"; TMT-11 and larger split it into 131N and 131C.
    """
    order = list(TMT_CHANNEL_MAP.keys())
    with open(out_path, "w") as f:
        for channel_key in sorted(channel_map.keys(),
                                  key=lambda c: order.index(c)
                                  if c in TMT_CHANNEL_MAP else 999):
            info = channel_map[channel_key]
            fp_channel = TMT_CHANNEL_MAP.get(channel_key, channel_key)
            if tmt_type == "TMT-10" and channel_key in ("tmt_131", "tmt_131n"):
                fp_channel = "131"
            sample = re.sub(r"\s+", "_", info["sample"].strip())
            f.write(f"{fp_channel}\t{sample}\n")


def patch_workflow(template_path, fasta_path, out_path, tmt_type):
    """
    Copy workflow template with per-plex fields patched for proteogenomics.

    Patches applied:
      - database.db-path: per-plex custom FASTA path
      - tmtintegrator.channel_num: detects and changes TMT
    """
    with open(template_path) as f:
        content = f.read()

    # Each entry: (key_regex, full_replacement_line, fallback_value_if_missing)
    patches = [
        (
            r"^database\.db-path=.*$",
            f"database.db-path={fasta_path}",
            f"database.db-path={fasta_path}",
        ),
        (
            r"^tmtintegrator\.channel_num=.*$",
            f"tmtintegrator.channel_num={tmt_type}",
            f"tmtintegrator.channel_num={tmt_type}",
        ),
    ]

    for key_pattern, replacement, fallback in patches:
        new_content = re.sub(key_pattern, replacement, content, flags=re.MULTILINE)
        if new_content == content:
            # Key not present in template — append it
            new_content = content.rstrip("\n") + f"\n{fallback}\n"
        content = new_content

    with open(out_path, "w") as f:
        f.write(content)


def main():
    parser = argparse.ArgumentParser(
        description="Generate FragPipe manifests for per-plex MS searches."
    )
    parser.add_argument("--tmt-map", required=True,
                        help="TMT mapping file (pdc_file_tmt_map.tsv)")
    parser.add_argument("--raw-dir", required=True,
                        help="Directory containing per-plex subdirectories "
                             "of RAW files (RAW_DIR/<plex>/*.raw)")
    parser.add_argument("--fasta-dir", default=None,
                        help="Directory containing per-plex FASTAs "
                             "(per-plex mode; mutually exclusive with "
                             "--single-fasta)")
    parser.add_argument("--single-fasta", default=None,
                        help="Path to one FASTA to search ALL plexes "
                             "against instead of per-plex FASTAs "
                             "(single mode; mutually exclusive with "
                             "--fasta-dir)")
    parser.add_argument("--workflow-template", default=None,
                        help="FragPipe workflow template to patch per plex")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for manifests and configs")

    args = parser.parse_args()

    if bool(args.fasta_dir) == bool(args.single_fasta):
        sys.exit("Error: specify exactly one of --fasta-dir (per-plex "
                 "mode) or --single-fasta (single mode)")

    for path, label in [(args.tmt_map, "TMT map"),
                        (args.raw_dir, "RAW directory")]:
        if not os.path.exists(path):
            sys.exit(f"Error: {label} not found: {path}")

    if args.fasta_dir and not os.path.exists(args.fasta_dir):
        sys.exit(f"Error: FASTA directory not found: {args.fasta_dir}")
    if args.single_fasta and not os.path.isfile(args.single_fasta):
        sys.exit(f"Error: Single FASTA not found: {args.single_fasta}")

    if args.workflow_template and not os.path.isfile(args.workflow_template):
        sys.exit(f"Error: Workflow template not found: "
                 f"{args.workflow_template}")

    # Create output directories
    manifest_dir = os.path.join(args.output, "manifests")
    workflow_dir = os.path.join(args.output, "workflows")
    results_dir = os.path.join(args.output, "results")
    for d in [manifest_dir, results_dir]:
        os.makedirs(d, exist_ok=True)
    if args.workflow_template:
        os.makedirs(workflow_dir, exist_ok=True)

    # Load TMT mapping
    print(f"Loading TMT map: {args.tmt_map}")
    plex_to_files, plex_to_channels = load_tmt_map(args.tmt_map)
    print(f"  {len(plex_to_files)} plexes")

    # Scan for available RAW files one level down (RAW_DIR/<plex>/<file>.raw)
    # Maps file name -> list of full paths (a name could appear in more
    # than one plex subdirectory)
    available_raws = defaultdict(list)
    n_subdirs = 0
    with os.scandir(args.raw_dir) as top:
        for subdir in top:
            if not subdir.is_dir():
                continue
            n_subdirs += 1
            with os.scandir(subdir.path) as inner:
                for entry in inner:
                    if entry.is_file() and entry.name.lower().endswith(".raw"):
                        available_raws[entry.name].append(entry.path)
    n_raw_paths = sum(len(p) for p in available_raws.values())
    print(f"  {n_raw_paths} RAW files in {n_subdirs} subdirectories "
          f"of {args.raw_dir}")
    duplicate_names = [n for n, p in available_raws.items() if len(p) > 1]
    if duplicate_names:
        print(f"  WARNING: {len(duplicate_names)} RAW file name(s) appear "
              f"in more than one subdirectory")

    # Scan for available per-plex FASTAs (per-plex mode only)
    available_fastas = set()
    if args.fasta_dir:
        for fname in os.listdir(args.fasta_dir):
            if fname.endswith(".fasta"):
                plex_id = fname.replace(".fasta", "")
                available_fastas.add(plex_id)
        print(f"  {len(available_fastas)} per-plex FASTAs available")
    else:
        print(f"  Single mode: all plexes will use {args.single_fasta}")

    # Process each plex
    print(f"\nGenerating manifests...")

    plex_list = []
    summary_rows = []
    total_skipped_no_fasta = 0
    total_skipped_no_raws = 0

    for plex_id in sorted(plex_to_files.keys()):
        expected_files = plex_to_files[plex_id]

        # Determine FASTA for this plex
        if args.single_fasta:
            fasta_path = args.single_fasta
        else:
            fasta_path = os.path.join(args.fasta_dir, f"{plex_id}.fasta")
            if plex_id not in available_fastas:
                total_skipped_no_fasta += 1
                continue

        # Find available RAW files for this plex
        found_raws = []
        missing_raws = []
        for fname in sorted(expected_files):
            paths = available_raws.get(fname)
            if not paths:
                missing_raws.append(fname)
                continue
            if len(paths) == 1:
                raw_path = paths[0]
            else:
                # Prefer the copy in a subdirectory named after this plex
                matching = [p for p in paths
                            if os.path.basename(os.path.dirname(p)) == plex_id]
                raw_path = matching[0] if matching else sorted(paths)[0]
                print(f"  WARNING: {fname} found in {len(paths)} "
                      f"subdirectories; using {raw_path}")
            found_raws.append(raw_path)

        if not found_raws:
            total_skipped_no_raws += 1
            continue

        # Write manifest
        manifest_path = os.path.join(manifest_dir,
                                     f"{plex_id}.fp-manifest")
        write_manifest(found_raws, plex_id, manifest_path)

        # Write annotation into the plex's RAW subdirectory, i.e. the folder
        # holding this plex's RAW files. Normally that is RAW_DIR/<plex_id>;
        # if the files are spread over several folders, prefer the one named
        # after the plex, otherwise the folder holding the most files.
        raw_parent_dirs = [os.path.dirname(p) for p in found_raws]
        unique_dirs = sorted(set(raw_parent_dirs))
        named_dir = os.path.join(args.raw_dir, plex_id)
        if len(unique_dirs) == 1:
            plex_raw_dir = unique_dirs[0]
        elif any(os.path.basename(d) == plex_id for d in unique_dirs):
            plex_raw_dir = next(d for d in unique_dirs
                                if os.path.basename(d) == plex_id)
        else:
            plex_raw_dir = max(unique_dirs, key=raw_parent_dirs.count)
        if len(unique_dirs) > 1:
            print(f"  WARNING: {plex_id} RAW files are in {len(unique_dirs)} "
                  f"folders; writing annotation to {plex_raw_dir}")
        elif os.path.normpath(plex_raw_dir) != os.path.normpath(named_dir):
            print(f"  NOTE: {plex_id} RAW files are in "
                  f"{os.path.basename(plex_raw_dir)}/, not {plex_id}/; "
                  f"writing annotation there")

        channel_map = get_unique_channel_annotation(
            plex_to_channels[plex_id], plex_id)
        tmt_type = detect_tmt_type(channel_map)
        annotation_path = os.path.join(plex_raw_dir,
                                       f"{plex_id}_annotation.txt")
        write_annotation(channel_map, plex_id, annotation_path, tmt_type)

        # Patch workflow if template provided
        if args.workflow_template:
            workflow_path = os.path.join(workflow_dir,
                                         f"{plex_id}.workflow")
            patch_workflow(args.workflow_template,
                          os.path.abspath(fasta_path),
                          workflow_path, tmt_type)

        # Create per-plex results directory
        os.makedirs(os.path.join(results_dir, plex_id), exist_ok=True)

        plex_list.append(plex_id)
        summary_rows.append({
            "plex_id": plex_id,
            "expected_raws": len(expected_files),
            "found_raws": len(found_raws),
            "missing_raws": len(missing_raws),
            "n_channels": len(channel_map),
            "tmt_type": tmt_type,
            "fasta": os.path.basename(fasta_path),
            "annotation": annotation_path,
        })

    # Write plex list for SLURM array
    plex_list_path = os.path.join(args.output, "plex_list.txt")
    with open(plex_list_path, "w") as f:
        for plex_id in plex_list:
            f.write(plex_id + "\n")

    # Write summary
    summary_path = os.path.join(args.output, "search_setup_summary.tsv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t",
                                fieldnames=["plex_id", "expected_raws",
                                            "found_raws", "missing_raws",
                                            "n_channels", "tmt_type", "fasta",
                                            "annotation"])
        writer.writeheader()
        writer.writerows(summary_rows)

    # Print summary
    total_raws = sum(r["found_raws"] for r in summary_rows)
    print(f"\n{'=' * 50}")
    print("MS Search Setup Summary")
    print(f"{'=' * 50}")
    print(f"Plexes ready:           {len(plex_list)}")
    print(f"Plexes skipped (no FASTA): {total_skipped_no_fasta}")
    print(f"Plexes skipped (no RAWs):  {total_skipped_no_raws}")
    print(f"Total RAW files:        {total_raws}")
    for t in ("TMT-10", "TMT-11"):
        n = sum(1 for r in summary_rows if r["tmt_type"] == t)
        print(f"  {t} plexes:           {n}")
    print(f"Plex list:              {plex_list_path}")
    print(f"Manifests:              {manifest_dir}")
    print(f"Annotations:            {args.raw_dir}/<plex>/<plex>_annotation.txt")
    if args.workflow_template:
        print(f"Workflows:              {workflow_dir}")
    print(f"Results dir:            {results_dir}")
    print(f"Summary:                {summary_path}")

    if any(r["missing_raws"] > 0 for r in summary_rows):
        print(f"\nWARNING: Some plexes have missing RAW files.")
        print("  These plexes will be searched with available fractions only.")
        for r in summary_rows:
            if r["missing_raws"] > 0:
                print(f"  {r['plex_id']}: {r['missing_raws']} missing "
                      f"of {r['expected_raws']} expected")


if __name__ == "__main__":
    main()
