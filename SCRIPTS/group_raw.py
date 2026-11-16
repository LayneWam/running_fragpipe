#!/usr/bin/env python3
"""
Group raw files in RAW_DIR into plex subdirectories inside RAW_DIR.

A file named   09CPTAC_BCprospective_Proteome_BC_20170117_f18.raw
is moved into  RAW_DIR/09CPTAC_BCprospective_Proteome_BC_20170117/

The plex name is the file name with the trailing fraction suffix
(_f<.>.raw, case-insensitive) removed.

Usage:
    python group_raw.py                 # move files
    python group_raw.py --dry-run       # show what would happen
    python group_raw.py --raw-dir /path/to/RAW
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

RAW_DIR = Path('/scratch/lam.wayn/RAW')

# <plex>_f<digits>.raw  ->  group(1) = plex
FRACTION_RE = re.compile(r'^(?P<plex>.+)_f.*\.raw$', re.IGNORECASE)


def plex_name(file_name):
    """Return the plex name for a raw file, or None if the name doesn't match."""
    m = FRACTION_RE.match(file_name)
    return m.group('plex') if m else None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--raw-dir', type=Path, default=RAW_DIR,
                        help=f'directory containing raw files (default: {RAW_DIR})')
    parser.add_argument('--dry-run', action='store_true',
                        help='print planned moves without moving anything')
    args = parser.parse_args()

    raw_dir = args.raw_dir
    if not raw_dir.is_dir():
        sys.exit(f'RAW directory not found: {raw_dir}')

    # Only top-level regular files; skip partial downloads (.raw.part etc.)
    raw_files = sorted(p for p in raw_dir.iterdir()
                       if p.is_file() and p.suffix.lower() == '.raw')

    by_plex = defaultdict(list)
    unmatched = []
    for p in raw_files:
        plex = plex_name(p.name)
        if plex:
            by_plex[plex].append(p)
        else:
            unmatched.append(p)

    moved = skipped = failed = 0
    for plex, files in sorted(by_plex.items()):
        plex_dir = raw_dir / plex
        print(f'{plex}: {len(files)} files')
        if not args.dry_run:
            plex_dir.mkdir(exist_ok=True)

        for src in files:
            dest = plex_dir / src.name
            if dest.exists():
                print(f'  [EXISTS] {dest} - not overwriting')
                skipped += 1
                continue
            if args.dry_run:
                print(f'  [DRY] {src.name} -> {plex}/')
                moved += 1
                continue
            try:
                src.rename(dest)  # same filesystem, so this is an instant move
                moved += 1
            except OSError as e:
                print(f'  [FAIL] {src}: {e}', file=sys.stderr)
                failed += 1

    if unmatched:
        print(f'\n{len(unmatched)} .raw files did not match <plex>_f<.>.raw and were left in place:')
        for p in unmatched:
            print(f'  {p.name}')

    action = 'Would move' if args.dry_run else 'Moved'
    print(f'\n{action}: {moved}   Already in place: {skipped}   Failed: {failed}   '
          f'Plexes: {len(by_plex)}')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
