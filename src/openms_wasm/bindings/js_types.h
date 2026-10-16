// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// JavaScript values with a TypeScript type: the plain objects and arrays the bindings
// return or take. They are emscripten::val, but declared as own types so that the
// generated openms.d.ts names their shape instead of 'any'. The definitions are
// registered in js_types.cpp (register_type).

#include <emscripten/val.h>

/// Like EMSCRIPTEN_DECLARE_VAL_TYPE, but implicitly constructible from val, so a binding
/// returns a plain val as its declared type.
#define OPENMS_WASM_JS_TYPE(name)                                   \
  struct name : public ::emscripten::val                            \
  {                                                                 \
    name(const ::emscripten::val& other) : ::emscripten::val(other) {} \
  };

namespace OpenMS::WASM::JS
{
  // generic
  OPENMS_WASM_JS_TYPE(NumberArray);     ///< number[] (results)
  OPENMS_WASM_JS_TYPE(NumberInput);     ///< ArrayLike<number> (arguments: arrays or typed arrays)
  OPENMS_WASM_JS_TYPE(StringArray);     ///< string[]
  OPENMS_WASM_JS_TYPE(MetaValue);       ///< a meta value (or undefined)
  OPENMS_WASM_JS_TYPE(MetaValueInput);  ///< a meta value to set
  OPENMS_WASM_JS_TYPE(MetaValues);      ///< {name: meta value}

  // parameters
  OPENMS_WASM_JS_TYPE(ParamValue);
  OPENMS_WASM_JS_TYPE(ParamValues);
  OPENMS_WASM_JS_TYPE(ParamEntry);
  OPENMS_WASM_JS_TYPE(ParamEntries);

  // spectra and chromatograms
  OPENMS_WASM_JS_TYPE(SpectrumPeaks);
  OPENMS_WASM_JS_TYPE(ChromatogramPeaks);
  OPENMS_WASM_JS_TYPE(PeakObject);
  OPENMS_WASM_JS_TYPE(OptionalPeak);
  OPENMS_WASM_JS_TYPE(BasePeak);
  OPENMS_WASM_JS_TYPE(PrecursorObject);
  OPENMS_WASM_JS_TYPE(PrecursorInput);
  OPENMS_WASM_JS_TYPE(PrecursorList);
  OPENMS_WASM_JS_TYPE(PrecursorInputList);
  OPENMS_WASM_JS_TYPE(ProductObject);
  OPENMS_WASM_JS_TYPE(ProductList);
  OPENMS_WASM_JS_TYPE(FloatDataArrays);
  OPENMS_WASM_JS_TYPE(IntegerDataArrays);
  OPENMS_WASM_JS_TYPE(StringDataArrays);
  OPENMS_WASM_JS_TYPE(SpectrumSummary);
  OPENMS_WASM_JS_TYPE(PeakData2D);

  // chemistry
  OPENMS_WASM_JS_TYPE(ResidueObject);
  OPENMS_WASM_JS_TYPE(ResidueList);
  OPENMS_WASM_JS_TYPE(ModificationObject);
  OPENMS_WASM_JS_TYPE(TermSpecificityInput);
  OPENMS_WASM_JS_TYPE(EnzymeObject);
  OPENMS_WASM_JS_TYPE(ElementCounts);
  OPENMS_WASM_JS_TYPE(IsotopePeaks);
  OPENMS_WASM_JS_TYPE(PeakObjectList);
  OPENMS_WASM_JS_TYPE(AASequenceList);

  // processing
  OPENMS_WASM_JS_TYPE(Alignment);
  OPENMS_WASM_JS_TYPE(ThreadingInfo);
  OPENMS_WASM_JS_TYPE(MemoryInfo);

  // files
  OPENMS_WASM_JS_TYPE(PeakFileOptions);
  OPENMS_WASM_JS_TYPE(FASTAEntries);

  // identifications
  OPENMS_WASM_JS_TYPE(PeptideHitObject);
  OPENMS_WASM_JS_TYPE(PeptideIdentificationObject);
  OPENMS_WASM_JS_TYPE(PeptideIdentificationObjects);
  OPENMS_WASM_JS_TYPE(ProteinHitObject);
  OPENMS_WASM_JS_TYPE(ProteinIdentificationObject);
  OPENMS_WASM_JS_TYPE(ProteinIdentificationObjects);
  OPENMS_WASM_JS_TYPE(PeptideHitList);
  OPENMS_WASM_JS_TYPE(ProteinHitList);
  OPENMS_WASM_JS_TYPE(PeakAnnotations);

  // ProSE
  OPENMS_WASM_JS_TYPE(ProgressCallback);
  OPENMS_WASM_JS_TYPE(ProSESummary);
  OPENMS_WASM_JS_TYPE(ProSEOptions);
  OPENMS_WASM_JS_TYPE(ProSEFilesOptions);
  OPENMS_WASM_JS_TYPE(ProSEFileSummary);
  OPENMS_WASM_JS_TYPE(ProSEMultiSummary);
  OPENMS_WASM_JS_TYPE(AnnotationOptions);
  OPENMS_WASM_JS_TYPE(SpectrumAnnotation);

  // features
  OPENMS_WASM_JS_TYPE(BoundingBox);
  OPENMS_WASM_JS_TYPE(OptionalBoundingBox);
  OPENMS_WASM_JS_TYPE(BoundingBoxList);
  OPENMS_WASM_JS_TYPE(HullPoints);
  OPENMS_WASM_JS_TYPE(FeatureObject);
  OPENMS_WASM_JS_TYPE(FeatureObjects);
  OPENMS_WASM_JS_TYPE(FeatureTable);
  OPENMS_WASM_JS_TYPE(ConsensusFeatureObject);
  OPENMS_WASM_JS_TYPE(ConsensusFeatureObjects);
  OPENMS_WASM_JS_TYPE(ColumnHeaders);
  OPENMS_WASM_JS_TYPE(MassTraceObjects);
} // namespace OpenMS::WASM::JS
