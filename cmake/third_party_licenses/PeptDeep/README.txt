PeptDeep models. share/OpenMS/models holds peptdeep_rt_dynamic.onnx, peptdeep_ms2_dynamic.onnx and
peptdeep_ccs_dynamic.onnx: the retention time, MS2 and CCS models of the pretrained generic models of
AlphaPeptDeep (https://github.com/MannLabs/alphapeptdeep), which OpenMS converted from PyTorch to
ONNX with tools/scripts/export_peptdeep_models_to_onnx.py. AlphaPeptDeep and its pretrained models
carry this notice:

  Copyright 2020 MannLabs

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

The Apache License 2.0 is in LICENSE.txt, next to this file. The conversion stores the trained
weights in another file format (torch.onnx.export, without quantization); the script that does it is
part of the OpenMS source code (https://github.com/OpenMS/OpenMS).
