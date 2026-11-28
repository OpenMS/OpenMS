// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Peptide and protein identifications (METADATA/) and the database search of
// SimpleSearchEngineAlgorithm (ANALYSIS/ID/).

#include "binding_utils.h"

#include <OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace emscripten;

namespace OpenMS::WASM
{
  using ProteinIdentificationList = std::vector<ProteinIdentification>;

  namespace
  {
    template <typename Container>
    void checkIndex(const Container& c, Size i)
    {
      if (i >= c.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, i, c.size());
    }

    //------------------------------------------------------------------------------
    // Plain objects (read-only views of the results)

    JS::PeptideHitObject peptideHitToJS(const PeptideHit& h)
    {
      val o = val::object();
      o.set("sequence", h.getSequence().toString());
      o.set("score", h.getScore());
      o.set("charge", h.getCharge());
      o.set("rank", h.getRank());
      o.set("isDecoy", h.isDecoy());
      o.set("proteinAccessions", toArray(h.extractProteinAccessionsSet()));
      o.set("metaValues", metaValuesToObject(h));
      return o;
    }

    JS::PeptideIdentificationObject peptideIdToJS(const PeptideIdentification& id)
    {
      val o = val::object();
      o.set("rt", id.hasRT() ? val(id.getRT()) : val::null());
      o.set("mz", id.hasMZ() ? val(id.getMZ()) : val::null());
      o.set("scoreType", id.getScoreType());
      o.set("higherScoreBetter", id.isHigherScoreBetter());
      o.set("identifier", id.getIdentifier());
      o.set("spectrumReference", id.getSpectrumReference());
      o.set("hits", toArray(id.getHits(), peptideHitToJS));
      o.set("metaValues", metaValuesToObject(id));
      return o;
    }

    JS::ProteinHitObject proteinHitToJS(const ProteinHit& h)
    {
      val o = val::object();
      o.set("accession", h.getAccession());
      o.set("score", h.getScore());
      o.set("description", h.getDescription());
      o.set("sequence", h.getSequence());
      o.set("coverage", h.getCoverage());
      o.set("metaValues", metaValuesToObject(h));
      return o;
    }

    JS::ProteinIdentificationObject proteinIdToJS(const ProteinIdentification& id)
    {
      val o = val::object();
      o.set("identifier", id.getIdentifier());
      o.set("searchEngine", id.getSearchEngine());
      o.set("searchEngineVersion", id.getSearchEngineVersion());
      o.set("scoreType", id.getScoreType());
      o.set("higherScoreBetter", id.isHigherScoreBetter());
      o.set("hits", toArray(id.getHits(), proteinHitToJS));
      o.set("metaValues", metaValuesToObject(id));
      return o;
    }

    //------------------------------------------------------------------------------
    // PeptideHit

    AASequence hitGetSequence(const PeptideHit& h) { return h.getSequence(); }
    void hitSetSequence(PeptideHit& h, const AASequence& s) { h.setSequence(s); }
    std::string hitGetSequenceString(const PeptideHit& h) { return h.getSequence().toString(); }
    void hitSetSequenceString(PeptideHit& h, const std::string& s) { h.setSequence(AASequence::fromString(s)); }
    void hitSetCharge(PeptideHit& h, int charge) { h.setCharge(charge); }
    JS::StringArray hitProteinAccessions(const PeptideHit& h) { return toArray(h.extractProteinAccessionsSet()); }
    /// The fragment ions matched to peaks of the (preprocessed) spectrum, as the search engine annotated them
    JS::PeakAnnotations hitPeakAnnotations(const PeptideHit& h)
    {
      return toArray(h.getPeakAnnotations(), [](const PeptideHit::PeakAnnotation& a) {
        val o = val::object();
        o.set("mz", a.mz);
        o.set("intensity", a.intensity);
        o.set("annotation", a.annotation);
        o.set("charge", a.charge);
        return o;
      });
    }

    //------------------------------------------------------------------------------
    // PeptideIdentification

    JS::PeptideHitList pepIdGetHits(const PeptideIdentification& id) { return toArray(id.getHits()); }
    void pepIdSetHits(PeptideIdentification& id, const JS::PeptideHitList& hits) { id.setHits(vectorFromJS<PeptideHit>(hits)); }
    void pepIdInsertHit(PeptideIdentification& id, const PeptideHit& h) { id.insertHit(h); }
    double pepIdHitCount(const PeptideIdentification& id) { return static_cast<double>(id.getHits().size()); }
    std::string pepIdScoreType(const PeptideIdentification& id) { return id.getScoreType(); }
    std::string pepIdIdentifier(const PeptideIdentification& id) { return id.getIdentifier(); }

    //------------------------------------------------------------------------------
    // PeptideIdentificationList

    double pepListSize(const PeptideIdentificationList& l) { return static_cast<double>(l.size()); }
    PeptideIdentification pepListGet(const PeptideIdentificationList& l, Size i) { checkIndex(l, i); return l[i]; }
    void pepListSet(PeptideIdentificationList& l, Size i, const PeptideIdentification& id) { checkIndex(l, i); l[i] = id; }
    void pepListPush(PeptideIdentificationList& l, const PeptideIdentification& id) { l.push_back(id); }
    void pepListClear(PeptideIdentificationList& l) { l.clear(); }
    JS::PeptideIdentificationObjects pepListToArray(const PeptideIdentificationList& l) { return toArray(l, peptideIdToJS); }

    //------------------------------------------------------------------------------
    // ProteinHit / ProteinIdentification

    std::string protHitAccession(const ProteinHit& h) { return h.getAccession(); }
    std::string protHitSequence(const ProteinHit& h) { return h.getSequence(); }
    void protHitSetSequence(ProteinHit& h, const std::string& s) { h.setSequence(s); }

    JS::ProteinHitList protIdGetHits(const ProteinIdentification& id) { return toArray(id.getHits()); }
    void protIdSetHits(ProteinIdentification& id, const JS::ProteinHitList& hits) { id.setHits(vectorFromJS<ProteinHit>(hits)); }
    void protIdInsertHit(ProteinIdentification& id, const ProteinHit& h) { id.insertHit(h); }
    double protIdHitCount(const ProteinIdentification& id) { return static_cast<double>(id.getHits().size()); }
    std::string protIdIdentifier(const ProteinIdentification& id) { return id.getIdentifier(); }
    std::string protIdSearchEngine(const ProteinIdentification& id) { return id.getSearchEngine(); }
    std::string protIdSearchEngineVersion(const ProteinIdentification& id) { return id.getSearchEngineVersion(); }
    std::string protIdScoreType(const ProteinIdentification& id) { return id.getScoreType(); }

    double protListSize(const ProteinIdentificationList& l) { return static_cast<double>(l.size()); }
    ProteinIdentification protListGet(const ProteinIdentificationList& l, Size i) { checkIndex(l, i); return l[i]; }
    void protListSet(ProteinIdentificationList& l, Size i, const ProteinIdentification& id) { checkIndex(l, i); l[i] = id; }
    void protListPush(ProteinIdentificationList& l, const ProteinIdentification& id) { l.push_back(id); }
    void protListClear(ProteinIdentificationList& l) { l.clear(); }
    JS::ProteinIdentificationObjects protListToArray(const ProteinIdentificationList& l) { return toArray(l, proteinIdToJS); }

    //------------------------------------------------------------------------------
    // SimpleSearchEngineAlgorithm

    std::string searchRun(const SimpleSearchEngineAlgorithm& s, const std::string& spectra, const std::string& database,
                          ProteinIdentificationList& proteins, PeptideIdentificationList& peptides)
    {
      switch (s.search(spectra, database, proteins, peptides))
      {
        case SimpleSearchEngineAlgorithm::ExitCodes::EXECUTION_OK: return "EXECUTION_OK";
        case SimpleSearchEngineAlgorithm::ExitCodes::INPUT_FILE_EMPTY: return "INPUT_FILE_EMPTY";
        case SimpleSearchEngineAlgorithm::ExitCodes::UNEXPECTED_RESULT: return "UNEXPECTED_RESULT";
        case SimpleSearchEngineAlgorithm::ExitCodes::ILLEGAL_PARAMETERS: return "ILLEGAL_PARAMETERS";
        default: return "UNKNOWN_ERROR";
      }
    }
  } // namespace
} // namespace OpenMS::WASM

using namespace OpenMS;
using namespace OpenMS::WASM;

EMSCRIPTEN_BINDINGS(openms_identification)
{
  class_<PeptideHit>("PeptideHit")
    .constructor<>()
    .constructor<const PeptideHit&>()
    .function("getSequence", OMS_FN(&hitGetSequence))
    .function("setSequence", OMS_FN(&hitSetSequence))
    .function("getSequenceString", OMS_FN(&hitGetSequenceString))
    .function("setSequenceString", OMS_FN(&hitSetSequenceString))
    .function("getScore", OMS_METHOD(PeptideHit, &PeptideHit::getScore))
    .function("setScore", OMS_METHOD(PeptideHit, &PeptideHit::setScore))
    .function("getCharge", OMS_METHOD(PeptideHit, &PeptideHit::getCharge))
    .function("setCharge", OMS_FN(&hitSetCharge))
    .function("getRank", OMS_METHOD(PeptideHit, &PeptideHit::getRank))
    .function("setRank", OMS_METHOD(PeptideHit, &PeptideHit::setRank))
    .function("isDecoy", OMS_METHOD(PeptideHit, &PeptideHit::isDecoy))
    .function("getProteinAccessions", OMS_FN(&hitProteinAccessions))
    .function("getPeakAnnotations", OMS_FN(&hitPeakAnnotations))
    .function("getMetaValue", OMS_FN(&metaValue<PeptideHit>))
    .function("setMetaValue", OMS_FN(&setMetaValue<PeptideHit>))
    .function("getMetaValues", OMS_FN(&metaValuesToObject<PeptideHit>))
    .function("toObject", OMS_FN(&peptideHitToJS))
    ;

  class_<PeptideIdentification>("PeptideIdentification")
    .constructor<>()
    .constructor<const PeptideIdentification&>()
    .function("getRT", OMS_METHOD(PeptideIdentification, &PeptideIdentification::getRT))
    .function("setRT", OMS_METHOD(PeptideIdentification, &PeptideIdentification::setRT))
    .function("hasRT", OMS_METHOD(PeptideIdentification, &PeptideIdentification::hasRT))
    .function("getMZ", OMS_METHOD(PeptideIdentification, &PeptideIdentification::getMZ))
    .function("setMZ", OMS_METHOD(PeptideIdentification, &PeptideIdentification::setMZ))
    .function("hasMZ", OMS_METHOD(PeptideIdentification, &PeptideIdentification::hasMZ))
    .function("getHits", OMS_FN(&pepIdGetHits))
    .function("setHits", OMS_FN(&pepIdSetHits))
    .function("insertHit", OMS_FN(&pepIdInsertHit))
    .function("getHitCount", OMS_FN(&pepIdHitCount))
    .function("getScoreType", OMS_FN(&pepIdScoreType))
    .function("setScoreType", OMS_METHOD(PeptideIdentification, &PeptideIdentification::setScoreType))
    .function("isHigherScoreBetter", OMS_METHOD(PeptideIdentification, &PeptideIdentification::isHigherScoreBetter))
    .function("setHigherScoreBetter", OMS_METHOD(PeptideIdentification, &PeptideIdentification::setHigherScoreBetter))
    .function("getIdentifier", OMS_FN(&pepIdIdentifier))
    .function("setIdentifier", OMS_METHOD(PeptideIdentification, &PeptideIdentification::setIdentifier))
    .function("getSpectrumReference", OMS_METHOD(PeptideIdentification, &PeptideIdentification::getSpectrumReference))
    .function("setSpectrumReference", OMS_METHOD(PeptideIdentification, &PeptideIdentification::setSpectrumReference))
    .function("sort", OMS_METHOD(PeptideIdentification, &PeptideIdentification::sort))
    .function("getMetaValue", OMS_FN(&metaValue<PeptideIdentification>))
    .function("setMetaValue", OMS_FN(&setMetaValue<PeptideIdentification>))
    .function("getMetaValues", OMS_FN(&metaValuesToObject<PeptideIdentification>))
    .function("toObject", OMS_FN(&peptideIdToJS))
    ;

  class_<PeptideIdentificationList>("PeptideIdentificationList")
    .constructor<>()
    .function("size", OMS_FN(&pepListSize))
    .function("get", OMS_FN(&pepListGet))
    .function("set", OMS_FN(&pepListSet))
    .function("push_back", OMS_FN(&pepListPush))
    .function("clear", OMS_FN(&pepListClear))
    .function("toArray", OMS_FN(&pepListToArray))
    ;

  class_<ProteinHit>("ProteinHit")
    .constructor<>()
    .constructor<const ProteinHit&>()
    .function("getAccession", OMS_FN(&protHitAccession))
    .function("setAccession", OMS_METHOD(ProteinHit, &ProteinHit::setAccession))
    .function("getScore", OMS_METHOD(ProteinHit, &ProteinHit::getScore))
    .function("setScore", OMS_METHOD(ProteinHit, &ProteinHit::setScore))
    .function("getSequence", OMS_FN(&protHitSequence))
    .function("setSequence", OMS_FN(&protHitSetSequence))
    .function("getDescription", OMS_METHOD(ProteinHit, &ProteinHit::getDescription))
    .function("setDescription", OMS_METHOD(ProteinHit, &ProteinHit::setDescription))
    .function("getCoverage", OMS_METHOD(ProteinHit, &ProteinHit::getCoverage))
    .function("setCoverage", OMS_METHOD(ProteinHit, &ProteinHit::setCoverage))
    .function("getMetaValue", OMS_FN(&metaValue<ProteinHit>))
    .function("setMetaValue", OMS_FN(&setMetaValue<ProteinHit>))
    .function("getMetaValues", OMS_FN(&metaValuesToObject<ProteinHit>))
    .function("toObject", OMS_FN(&proteinHitToJS))
    ;

  class_<ProteinIdentification>("ProteinIdentification")
    .constructor<>()
    .constructor<const ProteinIdentification&>()
    .function("getHits", OMS_FN(&protIdGetHits))
    .function("setHits", OMS_FN(&protIdSetHits))
    .function("insertHit", OMS_FN(&protIdInsertHit))
    .function("getHitCount", OMS_FN(&protIdHitCount))
    .function("getIdentifier", OMS_FN(&protIdIdentifier))
    .function("setIdentifier", OMS_METHOD(ProteinIdentification, &ProteinIdentification::setIdentifier))
    .function("getSearchEngine", OMS_FN(&protIdSearchEngine))
    .function("setSearchEngine", OMS_METHOD(ProteinIdentification, &ProteinIdentification::setSearchEngine))
    .function("getSearchEngineVersion", OMS_FN(&protIdSearchEngineVersion))
    .function("getScoreType", OMS_FN(&protIdScoreType))
    .function("setScoreType", OMS_METHOD(ProteinIdentification, &ProteinIdentification::setScoreType))
    .function("isHigherScoreBetter", OMS_METHOD(ProteinIdentification, &ProteinIdentification::isHigherScoreBetter))
    .function("setHigherScoreBetter", OMS_METHOD(ProteinIdentification, &ProteinIdentification::setHigherScoreBetter))
    .function("getMetaValue", OMS_FN(&metaValue<ProteinIdentification>))
    .function("setMetaValue", OMS_FN(&setMetaValue<ProteinIdentification>))
    .function("getMetaValues", OMS_FN(&metaValuesToObject<ProteinIdentification>))
    .function("toObject", OMS_FN(&proteinIdToJS))
    ;

  class_<ProteinIdentificationList>("ProteinIdentificationList")
    .constructor<>()
    .function("size", OMS_FN(&protListSize))
    .function("get", OMS_FN(&protListGet))
    .function("set", OMS_FN(&protListSet))
    .function("push_back", OMS_FN(&protListPush))
    .function("clear", OMS_FN(&protListClear))
    .function("toArray", OMS_FN(&protListToArray))
    ;

  class_<SimpleSearchEngineAlgorithm, base<DefaultParamHandler>>("SimpleSearchEngineAlgorithm")
    .constructor<>()
    .function("search", OMS_FN(&searchRun))
    ;
}
