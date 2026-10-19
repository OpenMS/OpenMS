#!/usr/bin/env python3
"""Build a 1:1 SHUFFLED entrapment database — the design behind the published glyco numbers.

    build_shuffled_entrap.py <targets.fasta> <out.fasta> [seed]

For every target sequence this appends a companion whose amino-acid composition and length
are identical but whose order is shuffled, tagged `ENTRAP_` inside the accession AND the
entry name (`>sp|ENTRAP_Q99JY4|ENTRAP_TRABD_MOUSE`), which is the layout `eval_entrap.py`
detects by substring.

WHY SHUFFLED-SELF RATHER THAN A FOREIGN PROTEOME
`build_entrap.py` in this directory does something different: it appends an unrelated
proteome (yeast / E. coli). Both are valid entrapment designs, but they are NOT
interchangeable for glyco, and mixing them is why an independent reproduction of the mouse
benchmark came out ~21% high:

  * Ratio. Shuffled-self is exactly 1:1, so the FDP correction factor is 2. A foreign
    proteome is whatever size it happens to be — mouse + E. coli reviewed is about 3.9:1,
    needing a factor near 4.9, and no one is reminded of that by the file name.
  * Sequon density. Only peptides carrying an N-X-S/T sequon can be glycopeptides, so the
    entrapment space that matters is sequon-bearing peptides, not all peptides. Shuffling
    preserves composition, so sequon density tracks the target closely. A foreign proteome
    does not, which forces a large and easily-forgotten sequon correction.
  * Search-space size. A smaller database yields fewer candidates and fewer decoys, which
    moves the q-value threshold and therefore the count — independently of any real change
    in engine behaviour.

Use this script to reproduce the published mouse glyco numbers. Use
`build_entrap.py` when you deliberately want foreign-proteome entrapment, and then size the
correction factor from the database you actually built rather than assuming 1:1.
"""
import random
import sys


def read_fasta(path):
    hdr, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(buf)
            hdr, buf = line.rstrip("\n"), []
        else:
            buf.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(buf)


def tag(hdr):
    """`>sp|Q99JY4|TRABD_MOUSE ...` -> `>sp|ENTRAP_Q99JY4|ENTRAP_TRABD_MOUSE ...`"""
    parts = hdr.split("|")
    if len(parts) >= 3:
        parts[1] = "ENTRAP_" + parts[1]
        parts[2] = "ENTRAP_" + parts[2]
        return "|".join(parts)
    return ">ENTRAP_" + hdr.lstrip(">")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    src, out = sys.argv[1], sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = random.Random(seed)  # deterministic: same input + seed -> same database

    n = 0
    with open(out, "w") as fh:
        entries = list(read_fasta(src))
        for hdr, seq in entries:                      # targets first, verbatim
            fh.write(f"{hdr}\n{seq}\n")
        for hdr, seq in entries:                      # then the shuffled twins
            chars = list(seq)
            rng.shuffle(chars)
            fh.write(f"{tag(hdr)}\n{''.join(chars)}\n")
            n += 1
    print(f"wrote {out}: {len(entries)} targets + {n} shuffled entrapment (1:1, seed={seed})")
    print("FDP correction factor for this database = 2.0 (it is exactly 1:1)")
    print("run andes with --decoy-strategy sequon-reverse so decoys are built correctly")


if __name__ == "__main__":
    main()
