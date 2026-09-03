// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief The label state of an identification in MS1-labeled (SILAC, Dimethyl, ...) data

    In MS1-labeled data the label belongs to the channel, not to the peptide: the light and the heavy
    spectrum of one peptide have to name one peptide for linking, match between runs, inference and
    quantification. MS1LabeledWorkflow therefore reduces the sequence of every identification to that
    <em>peptide identity</em> (the label modifications removed) and records the label state on the
    PeptideHit as meta values:

      - @c labeled_sequence: the peptidoform the spectrum was matched with, e.g.
        <tt>PEPTIDEK(Label:13C(6)15N(2))</tt>
      - @c removed_labels: the labels removed from it, in the FeatureFinderMultiplex vocabulary
        (e.g. @c Lys8, or @c Arg10,Lys8), or @c none
      - @c label_channel: the 1-based channel the spectrum belongs to, i.e. the @c Label of the
        experimental design (1 = light); 0 if it could not be determined

    PSM-level output describes the spectrum match and therefore reports the matched peptidoform
    (mzTab PSM section, QPX psm view), so that the calculated mass fits the precursor; feature-level
    output reports the peptide identity. The QPX feature and psm views additionally carry the three
    values as @c cv_params.

    @ingroup Metadata
  */
  namespace MS1LabelState
  {
    /// Meta value key: the peptidoform the spectrum was matched with (see the namespace documentation)
    inline const std::string LABELED_SEQUENCE = "labeled_sequence";
    /// Meta value key: the labels removed from the matched peptidoform, or "none"
    inline const std::string REMOVED_LABELS = "removed_labels";
    /// Meta value key: the 1-based channel the spectrum belongs to (0 = unknown)
    inline const std::string LABEL_CHANNEL = "label_channel";

    /// The three keys, in the order they are reported
    OPENMS_DLLAPI const std::vector<std::string>& keys();

    /**
      @brief The registry indices of the three keys, resolved once for a loop over many hits

      A meta value looked up by name takes the registry lock on every call; the exporters walk
      millions of hits, so they resolve the indices once per table and look up by index. A key that
      no hit in this process has used yet resolves to @c UInt(-1), which the index overloads treat
      as "absent". Construct it right before the loop, not statically: the keys become registered
      the moment a workflow annotates its first hit.
    */
    struct OPENMS_DLLAPI Keys
    {
      Keys();
      UInt labeled_sequence;
      UInt removed_labels;
      UInt label_channel;
    };

    /// Does @p hit record a non-empty @c labeled_sequence, i.e. was it reduced to a peptide identity?
    /// Presence only: the recorded peptidoform may equal the hit's own sequence, as it does for an
    /// unlabeled channel. An empty value counts as absent, since it records no peptidoform.
    OPENMS_DLLAPI bool hasMatchedSequence(const PeptideHit& hit);
    /// Index-based variant of hasMatchedSequence() for loops over many hits
    OPENMS_DLLAPI bool hasMatchedSequence(const PeptideHit& hit, const Keys& keys);

    /**
      @brief The peptidoform @p hit was matched to the spectrum with

      The hit's sequence, unless the hit records a @c labeled_sequence, in which case that one.

      Never throws. The value comes from a file and may be absent, empty or unparsable; any of those
      yields the hit's own sequence, with a warning logged once per process for an unparsable one.
      Callers stream their output, so failing a whole export on one bad row would abandon a
      half-written file.
    */
    OPENMS_DLLAPI AASequence matchedSequence(const PeptideHit& hit);
    /// Index-based variant of matchedSequence() for loops over many hits
    OPENMS_DLLAPI AASequence matchedSequence(const PeptideHit& hit, const Keys& keys);

    /// A copy of @p hit whose sequence is matchedSequence(); an unchanged copy when no label state is recorded
    OPENMS_DLLAPI PeptideHit withMatchedSequence(const PeptideHit& hit);
  } // namespace MS1LabelState
} // namespace OpenMS
