#!/usr/bin/env python3
"""
Group raw files in RAW_DIR into plex subdirectories, using the PDC metadata.

Each raw file is moved to RAW_DIR/<plex_id>/, from run_metadata_id
of the pdc_file_tmt_map

Raw files are found directly in RAW_DIR and one folder down, so files that
were grouped earlier under a different folder name are moved to the right
plex folder. Old folders left empty are removed; folders that still hold
anything else (e.g. old annotation files) are left alone and reported.

Usage:
    python group_raw.py --dry-run       # show what would happen
    python group_raw.py                 # move files
    python group_raw.py --raw-dir /path/to/RAW --tmt-map /path/to/pdc_file_tmt_map.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

RAW_DIR = Path('/scratch/lam.wayn/RAW')
TMT_MAP = Path('/home/lam.wayn/aas_stabilization/reprocessing/PDC_meta/pdc_file_tmt_map.tsv')


def load_file_to_plex(tmt_map):
    """Map file_name -> run_metadata_id from pdc_file_tmt_map.tsv."""
    file_to_plexes = defaultdict(set)
    with open(tmt_map, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        missing = {'file_name', 'run_metadata_id'} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f'{tmt_map} is missing column(s): {sorted(missing)}')
        for row in reader:
            name = row['file_name'].strip()
            plex = row['run_metadata_id'].strip()
            if name and plex:
                file_to_plexes[name].add(plex)

    conflicts = {n: p for n, p in file_to_plexes.items() if len(p) > 1}
    if conflicts:
        print(f'WARNING: {len(conflicts)} file(s) belong to more than one plex '
              f'in {tmt_map.name}; they will be left in place:')
        for n, p in sorted(conflicts.items()):
            print(f'  {n}: {sorted(p)}')

    return {n: next(iter(p)) for n, p in file_to_plexes.items() if len(p) == 1}


def find_raw_files(raw_dir):
    """Raw files directly in raw_dir and one folder down (partial downloads excluded)."""
    found = []
    for entry in raw_dir.iterdir():
        if entry.is_file() and entry.suffix.lower() == '.raw':
            found.append(entry)
        elif entry.is_dir():
            found.extend(p for p in entry.iterdir()
                         if p.is_file() and p.suffix.lower() == '.raw')
    return sorted(found)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--raw-dir', type=Path, default=RAW_DIR,
                        help=f'directory containing raw files (default: {RAW_DIR})')
    parser.add_argument('--tmt-map', type=Path, default=TMT_MAP,
                        help=f'pdc_file_tmt_map.tsv (default: {TMT_MAP})')
    parser.add_argument('--dry-run', action='store_true',
                        help='print planned moves without moving anything')
    args = parser.parse_args()

    raw_dir = args.raw_dir
    if not raw_dir.is_dir():
        sys.exit(f'RAW directory not found: {raw_dir}')
    if not args.tmt_map.is_file():
        sys.exit(f'TMT map not found: {args.tmt_map}')

    file_to_plex = load_file_to_plex(args.tmt_map)
    print(f'{len(file_to_plex)} files assigned to '
          f'{len(set(file_to_plex.values()))} plexes in {args.tmt_map.name}\n')

    by_plex = defaultdict(list)
    unmapped = []
    for p in find_raw_files(raw_dir):
        plex = file_to_plex.get(p.name)
        if plex:
            by_plex[plex].append(p)
        else:
            unmapped.append(p)

    moved = in_place = conflicts = failed = 0
    old_dirs = set()

    for plex, files in sorted(by_plex.items()):
        safe_plex = plex.replace('/', '_')
        plex_dir = raw_dir / safe_plex
        to_move = [p for p in files if p.parent != plex_dir]
        in_place += len(files) - len(to_move)
        print(f'{safe_plex}: {len(files)} files ({len(to_move)} to move)')
        if not to_move:
            continue
        if not args.dry_run:
            plex_dir.mkdir(exist_ok=True)

        for src in to_move:
            dest = plex_dir / src.name
            if dest.exists():
                print(f'  [EXISTS] {dest} already present; leaving {src} in place')
                conflicts += 1
                continue
            where = src.parent.name + '/' if src.parent != raw_dir else 'RAW/'
            if args.dry_run:
                print(f'  [DRY] {where}{src.name} -> {safe_plex}/')
                moved += 1
                continue
            try:
                src.rename(dest)  # same filesystem: instant move, no copy
                moved += 1
                if src.parent != raw_dir:
                    old_dirs.add(src.parent)
            except OSError as e:
                print(f'  [FAIL] {src}: {e}', file=sys.stderr)
                failed += 1

    # Tidy up folders emptied by the moves
    leftover = []
    for d in sorted(old_dirs):
        try:
            d.rmdir()  # only succeeds if empty
        except OSError:
            leftover.append(d)

    if unmapped:
        print(f'\n{len(unmapped)} raw file(s) are not in {args.tmt_map.name} and were left in place:')
        for p in unmapped:
            print(f'  {p.relative_to(raw_dir)}')

    if leftover:
        print(f'\n{len(leftover)} old folder(s) still contain other files (e.g. old annotations):')
        for d in leftover:
            print(f'  {d.name}/: {", ".join(sorted(x.name for x in d.iterdir()))}')

    action = 'Would move' if args.dry_run else 'Moved'
    print(f'\n{action}: {moved}   Already in plex folder: {in_place}   '
          f'Name clashes: {conflicts}   Failed: {failed}   Plexes: {len(by_plex)}')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
