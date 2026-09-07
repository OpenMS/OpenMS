#!/usr/bin/env python3

import argparse
import importlib
from pathlib import Path

import numpy as np
import onnx
import onnxruntime as ort
import torch

from peptdeep.model.rt import AlphaRTModel
from peptdeep.settings import model_const


DEFAULT_PRETRAINED_DIR = Path(
    "~/peptdeep/pretrained_models/pretrained_models_v3/generic"
)
DEFAULT_OUT_DIR = Path(".")

MOD_HIDDEN = len(model_const["mod_elements"])


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export PeptDeep RT, MS2, and CCS pretrained PyTorch models to ONNX."
    )
    parser.add_argument(
        "--pretrained-dir",
        type=Path,
        default=DEFAULT_PRETRAINED_DIR,
        help=f"Directory containing rt.pth, ms2.pth, and ccs.pth. Default: {DEFAULT_PRETRAINED_DIR}",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help=f"Directory to write ONNX files. Default: {DEFAULT_OUT_DIR}",
    )
    return parser.parse_args()


def get_model_class(module_name, candidate_names):
    module = importlib.import_module(module_name)
    for name in candidate_names:
        if hasattr(module, name):
            return getattr(module, name)

    available = [name for name in dir(module) if "Model" in name or "model" in name]
    raise ImportError(
        f"Could not find any of {candidate_names} in {module_name}. "
        f"Available model-like names: {available}"
    )


def make_dummy_inputs(task, batch_size=1, seq_len=16):
    aa = torch.zeros((batch_size, seq_len), dtype=torch.long)
    mod_x = torch.zeros((batch_size, seq_len, MOD_HIDDEN), dtype=torch.float32)
    charges = torch.full((batch_size, 1), 2.0, dtype=torch.float32)
    nce = torch.full((batch_size, 1), 30.0, dtype=torch.float32)
    instrument_indices = torch.zeros((batch_size,), dtype=torch.long)

    if task == "rt":
        return aa, mod_x

    if task == "ccs":
        return aa, mod_x, charges

    if task == "ms2":
        return aa, mod_x, charges, nce, instrument_indices

    raise ValueError(f"Unknown PeptDeep task: {task}")


def warmup_and_load(model, pth_path, task):
    model.set_device("cpu")
    model.model.to(model.device)
    model.model.eval()

    with torch.no_grad():
        model.model(*make_dummy_inputs(task, seq_len=6))

    model.load(str(pth_path))
    model.model.eval()
    return model


def export_and_check(model, args, out_path, input_names, output_names, dynamic_axes):
    out_path = Path(out_path)

    torch.onnx.export(
        model.model,
        args,
        str(out_path),
        input_names=input_names,
        output_names=output_names,
        dynamic_axes=dynamic_axes,
        opset_version=17,
        do_constant_folding=True,
        dynamo=False,
        external_data=False,
    )

    onnx_model = onnx.load(str(out_path))
    onnx.checker.check_model(onnx_model)

    print(f"\nWrote {out_path}")
    sess = ort.InferenceSession(str(out_path), providers=["CPUExecutionProvider"])

    print("Inputs:")
    for inp in sess.get_inputs():
        print(" ", inp.name, inp.type, inp.shape)

    print("Outputs:")
    for out in sess.get_outputs():
        print(" ", out.name, out.type, out.shape)

    return sess


def export_rt(pretrained_dir, out_dir):
    model = AlphaRTModel()
    model = warmup_and_load(model, pretrained_dir / "rt.pth", "rt")

    return export_and_check(
        model=model,
        args=make_dummy_inputs("rt"),
        out_path=out_dir / "peptdeep_rt_dynamic.onnx",
        input_names=["input_sequences", "mod_x"],
        output_names=["rt_pred"],
        dynamic_axes={
            "input_sequences": {0: "batch_size", 1: "seq_length"},
            "mod_x": {0: "batch_size", 1: "seq_length"},
            "rt_pred": {0: "batch_size"},
        },
    )


def export_ms2(pretrained_dir, out_dir):
    AlphaMS2Model = get_model_class(
        "peptdeep.model.ms2",
        ["AlphaMS2Model", "pDeepModel", "MS2Model"],
    )

    model = AlphaMS2Model()
    model = warmup_and_load(model, pretrained_dir / "ms2.pth", "ms2")

    return export_and_check(
        model=model,
        args=make_dummy_inputs("ms2"),
        out_path=out_dir / "peptdeep_ms2_dynamic.onnx",
        input_names=["aa_indices", "mod_x", "charges", "nce", "instrument_indices"],
        output_names=["ms2_intensities"],
        dynamic_axes={
            "aa_indices": {0: "batch_size", 1: "seq_length"},
            "mod_x": {0: "batch_size", 1: "seq_length"},
            "charges": {0: "batch_size"},
            "nce": {0: "batch_size"},
            "instrument_indices": {0: "batch_size"},
            "ms2_intensities": {0: "batch_size", 1: "frag_position"},
        },
    )


def export_ccs(pretrained_dir, out_dir):
    AlphaCCSModel = get_model_class(
        "peptdeep.model.ccs",
        ["AlphaCCSModel", "CCSModel"],
    )

    model = AlphaCCSModel()
    model = warmup_and_load(model, pretrained_dir / "ccs.pth", "ccs")

    return export_and_check(
        model=model,
        args=make_dummy_inputs("ccs"),
        out_path=out_dir / "peptdeep_ccs_dynamic.onnx",
        input_names=["aa_indices", "mod_x", "charges"],
        output_names=["ccs_pred"],
        dynamic_axes={
            "aa_indices": {0: "batch_size", 1: "seq_length"},
            "mod_x": {0: "batch_size", 1: "seq_length"},
            "charges": {0: "batch_size"},
            "ccs_pred": {0: "batch_size"},
        },
    )


def smoke_rt(sess):
    aa = np.zeros((2, 11), dtype=np.int64)
    mod_x = np.zeros((2, 11, MOD_HIDDEN), dtype=np.float32)

    y = sess.run(None, {"input_sequences": aa, "mod_x": mod_x})[0]
    print("\nRT smoke:", y.shape, y.reshape(-1)[:3])


def smoke_ms2(sess):
    aa = np.zeros((2, 11), dtype=np.int64)
    mod_x = np.zeros((2, 11, MOD_HIDDEN), dtype=np.float32)
    charges = np.array([[2.0], [3.0]], dtype=np.float32)
    nce = np.array([[30.0], [30.0]], dtype=np.float32)
    instrument_indices = np.array([0, 0], dtype=np.int64)

    y = sess.run(
        None,
        {
            "aa_indices": aa,
            "mod_x": mod_x,
            "charges": charges,
            "nce": nce,
            "instrument_indices": instrument_indices,
        },
    )[0]

    print("\nMS2 smoke:", y.shape)


def smoke_ccs(sess):
    aa = np.zeros((2, 11), dtype=np.int64)
    mod_x = np.zeros((2, 11, MOD_HIDDEN), dtype=np.float32)
    charges = np.array([[2.0], [3.0]], dtype=np.float32)

    y = sess.run(
        None,
        {
            "aa_indices": aa,
            "mod_x": mod_x,
            "charges": charges,
        },
    )[0]

    print("\nCCS smoke:", y.shape, y.reshape(-1)[:3])


if __name__ == "__main__":
    args = parse_args()

    torch.manual_seed(1337)

    pretrained_dir = args.pretrained_dir.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()

    if not pretrained_dir.exists():
        raise FileNotFoundError(
            f"Pretrained directory does not exist: {pretrained_dir}"
        )

    for filename in ["rt.pth", "ms2.pth", "ccs.pth"]:
        path = pretrained_dir / filename
        if not path.exists():
            raise FileNotFoundError(f"Missing pretrained model: {path}")

    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Using pretrained models from: {pretrained_dir}")
    print(f"Writing ONNX models to: {out_dir}")

    rt_sess = export_rt(pretrained_dir, out_dir)
    smoke_rt(rt_sess)

    ms2_sess = export_ms2(pretrained_dir, out_dir)
    smoke_ms2(ms2_sess)

    ccs_sess = export_ccs(pretrained_dir, out_dir)
    smoke_ccs(ccs_sess)
