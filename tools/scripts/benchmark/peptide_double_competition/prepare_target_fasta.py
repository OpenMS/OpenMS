#!/usr/bin/env python3
"""Keep only target proteins from the paired OpenMS tutorial FASTA."""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    count = 0
    write = False
    with args.input.open(encoding="ascii") as source, args.output.open("w", encoding="ascii", newline="\n") as target:
        for line in source:
            if line.startswith(">"):
                write = "DECOY_" not in line
                count += write
            if write:
                target.write(line)
    print(f"target_proteins={count} output={args.output}")


if __name__ == "__main__":
    main()
