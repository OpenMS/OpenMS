#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

"""Developer-only prediction parity check. Native analysis never invokes Python."""

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from train_trace_model import TraceCNN


def check(native_test, directory):
    subprocess.run([str(native_test), "--parity", str(directory)], check=True)
    # C++ OutputArchive stores a readable module, including its parameter names.
    archive = torch.jit.load(str(directory / "native.pt"), map_location="cpu")
    weights = {name.removeprefix("model."): value for name, value in archive.state_dict().items()
               if name.startswith("model.")}
    model = TraceCNN().eval()
    model.load_state_dict(weights, strict=True)
    assert sum(p.numel() for p in model.parameters()) == 6141
    raw = np.fromfile(directory / "raw.f32", dtype="<f4").reshape(-1, 6, 32, 8)
    auxiliary = np.fromfile(directory / "aux.f32", dtype="<f4").reshape(-1, 7)
    expected = np.fromfile(directory / "logits.f32", dtype="<f4")
    x = np.log1p(raw)
    x /= np.maximum(x.max(axis=(1, 2, 3), keepdims=True), 1e-12)
    with torch.no_grad():
        actual = model(torch.from_numpy(x), torch.from_numpy(auxiliary)).numpy()
    np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-5)
    print(f"Native/Python identical-weight prediction parity: {len(actual)} windows, "
          f"max_abs_difference={np.max(np.abs(actual - expected)):.9g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-test", required=True, type=Path)
    args = parser.parse_args()
    torch.set_num_threads(1)
    with tempfile.TemporaryDirectory(prefix="opendia-native-parity-") as temporary:
        check(args.native_test.resolve(), Path(temporary))
