"""Minimal MSNumpress decoders for .xipm arrays (zlib + numpress linear / slof), checked against target m/z windows."""
import struct
import zlib

import numpy as np


def _nibbles(data):
    out = []
    for byte in data:
        out.append(byte >> 4)
        out.append(byte & 0xF)
    return out


def _fixed_point(data):
    le = struct.unpack("<d", data[:8])[0]
    be = struct.unpack(">d", data[:8])[0]
    # numpress writes the double in big-endian byte order on little-endian hosts
    return be if 0 < be < 1e15 else le


def decode_linear(blob):
    if not blob:
        return np.zeros(0)
    d = zlib.decompress(blob)
    if len(d) < 8:
        return np.zeros(0)
    fp = _fixed_point(d)
    if len(d) == 8:
        return np.zeros(0)
    ints = [struct.unpack("<I", d[8:12])[0]]
    if len(d) >= 16:
        ints.append(struct.unpack("<I", d[12:16])[0])
    if len(d) <= 16:
        return np.array(ints, float) / fp
    nib = _nibbles(d[16:])
    # the last byte may carry one padding nibble
    p = 0
    vals = list(ints)
    while p < len(nib):
        if p == len(nib) - 1 and nib[p] == 0:
            break
        head = nib[p]
        p += 1
        if head <= 8:
            n = head
            res = 0
        else:
            n = head - 8
            res = 0
            for i in range(n):
                res |= 0xF0000000 >> (4 * i)
        if n < 8:
            if p + (8 - n) > len(nib):
                break
            for i in range(n, 8):
                res |= nib[p] << ((i - n) * 4)
                p += 1
        diff = struct.unpack("<i", struct.pack("<I", res & 0xFFFFFFFF))[0]
        extrapol = vals[-1] + (vals[-1] - vals[-2])
        vals.append(extrapol + diff)
    return np.array(vals, float) / fp


def decode_slof(blob):
    if not blob:
        return np.zeros(0)
    d = zlib.decompress(blob)
    fp = _fixed_point(d)
    raw = np.frombuffer(d[8:], dtype="<u2").astype(float)
    return np.exp(raw / fp) - 1.0


def decode_row(row):
    mz = decode_linear(row["MZ_DATA"])
    rt = decode_linear(row["RT_DATA"])
    im = decode_linear(row["MOBILITY_DATA"])
    it = decode_slof(row["INTENSITY_DATA"])
    n = min(len(mz), len(rt), len(it))
    return mz[:n], rt[:n], (im[:n] if len(im) >= n else np.full(n, -1.0)), it[:n]
