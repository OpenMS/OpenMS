// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MS1LabelState.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/MetaInfo.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>

#include <atomic>

namespace OpenMS::MS1LabelState
{
  const std::vector<std::string>& keys()
  {
    static const std::vector<std::string> k{LABELED_SEQUENCE, REMOVED_LABELS, LABEL_CHANNEL};
    return k;
  }

  Keys::Keys() :
    labeled_sequence(MetaInfo::registry().getIndex(LABELED_SEQUENCE)),
    removed_labels(MetaInfo::registry().getIndex(REMOVED_LABELS)),
    label_channel(MetaInfo::registry().getIndex(LABEL_CHANNEL))
  {
  }

  namespace
  {
    /**
      @brief The recorded peptidoform, or @p fallback when it cannot be used

      The value reaches us from a file: it survives consensusXML and idXML as a plain UserParam, so
      it can be edited, produced by a foreign tool, or left behind by one that rewrote the sequence.
      A hit that cannot be read must cost that hit its label state, not the whole export -- the
      exporters stream, so a throw here abandons a half-written file (and the empty-string case did
      not even throw: it silently replaced the sequence with an empty one and dropped the row).

      Exception::BaseException is caught rather than the individual types because
      AASequence::fromString has three reachable ones, and the least obvious is the likeliest here:
      square-bracket notation such as PEPTIDEK[UNIMOD:259], the form the QPX peptidoform column
      emits, raises ConversionError rather than ParseError or InvalidValue.

      Falling back to the identity can make the row collide with a sibling that legitimately carries
      it -- the light and heavy match of one spectrum share it. The QPX psm identity check then
      refuses the table, which is the right answer: two rows really do claim the same primary key.
    */
    AASequence parseOrFallback(const std::string& value, const AASequence& fallback)
    {
      try
      {
        return AASequence::fromString(value);
      }
      catch (const Exception::BaseException& e)
      {
        // Once per process: a broken file breaks every row of it, and two of the three call sites
        // run inside an OpenMP region (OPENMS_LOG_WARN itself is thread-safe, see LogStream.h).
        static std::atomic<bool> warned{false};
        if (!warned.exchange(true))
        {
          OPENMS_LOG_WARN << "MS1LabelState: the '" << LABELED_SEQUENCE << "' meta value '" << value
                          << "' is not a peptide sequence (" << e.getName() << ": " << e.getMessage()
                          << "). Reporting the hit's own sequence instead; further such values are "
                             "handled the same way without another warning." << std::endl;
        }
        return fallback;
      }
    }
  } // namespace

  bool hasMatchedSequence(const PeptideHit& hit)
  {
    // An empty value records no peptidoform: parsing it yields an empty AASequence rather than
    // throwing, so without this it would silently replace the hit's sequence with nothing.
    return hit.metaValueExists(LABELED_SEQUENCE) && !hit.getMetaValue(LABELED_SEQUENCE).toString().empty();
  }

  bool hasMatchedSequence(const PeptideHit& hit, const Keys& keys)
  {
    // UInt(-1) is guarded explicitly: the index overload of metaValueExists() does not special-case
    // the "not registered" sentinel the way the string overload does.
    return keys.labeled_sequence != static_cast<UInt>(-1) && hit.metaValueExists(keys.labeled_sequence)
           && !hit.getMetaValue(keys.labeled_sequence).toString().empty();
  }

  AASequence matchedSequence(const PeptideHit& hit, const Keys& keys)
  {
    if (!hasMatchedSequence(hit, keys)) { return hit.getSequence(); }
    return parseOrFallback(hit.getMetaValue(keys.labeled_sequence).toString(), hit.getSequence());
  }

  AASequence matchedSequence(const PeptideHit& hit)
  {
    if (!hasMatchedSequence(hit)) { return hit.getSequence(); }
    return parseOrFallback(hit.getMetaValue(LABELED_SEQUENCE).toString(), hit.getSequence());
  }

  PeptideHit withMatchedSequence(const PeptideHit& hit)
  {
    PeptideHit matched = hit;
    matched.setSequence(matchedSequence(hit));
    return matched;
  }
} // namespace OpenMS::MS1LabelState
