#!/usr/bin/env python3
"""Fetch a public file with parallel curl byte ranges and verify the result size."""

import argparse
import concurrent.futures
import subprocess
from pathlib import Path


def download(url, output, start, end):
    part = output.with_name(f"{output.name}.part{start}")
    expected = end - start + 1
    subprocess.run(
        ["curl.exe", "--fail", "--location", "--silent", "--show-error", "--retry", "3",
         "--range", f"{start}-{end}", "--output", str(part), url],
        check=True,
    )
    if part.stat().st_size != expected:
        raise RuntimeError(f"Wrong byte-range size for {start}-{end}: {part.stat().st_size} != {expected}")
    return part


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("url")
    parser.add_argument("output", type=Path)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    chunk_size = (args.size + args.workers - 1) // args.workers
    ranges = [(start, min(args.size - 1, start + chunk_size - 1)) for start in range(0, args.size, chunk_size)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(download, args.url, args.output, start, end) for start, end in ranges]
        parts = [future.result() for future in futures]
    with args.output.open("wb") as target:
        for part in parts:
            with part.open("rb") as source:
                while block := source.read(8 * 1024 * 1024):
                    target.write(block)
    if args.output.stat().st_size != args.size:
        raise RuntimeError(f"Wrong assembled file size: {args.output.stat().st_size} != {args.size}")
    for part in parts:
        part.unlink()
    print(f"Downloaded {args.output} ({args.size:,} bytes)")


if __name__ == "__main__":
    main()
