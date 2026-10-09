// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/METADATA/AnnotatedMSRun.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <functional>
#include <optional>
#include <unordered_set>

using namespace std;

namespace OpenMS
{

  IDMapper::IDMapper() :
    DefaultParamHandler("IDMapper"),
    rt_tolerance_(5.0),
    mz_tolerance_(20),
    measure_(MEASURE_PPM),
    ignore_charge_(false)
  {
    defaults_.setValue("rt_tolerance", rt_tolerance_, "RT tolerance (in seconds) for the matching");
    defaults_.setMinFloat("rt_tolerance", 0);
    defaults_.setValue("mz_tolerance", mz_tolerance_, "m/z tolerance (in ppm or Da) for the matching");
    defaults_.setMinFloat("mz_tolerance", 0);
    defaults_.setValue("mz_measure", "ppm", "unit of 'mz_tolerance' (ppm or Da)");
    defaults_.setValidStrings("mz_measure", {"ppm","Da"});
    defaults_.setValue("mz_reference", "precursor", "source of m/z values for peptide identifications");
    defaults_.setValidStrings("mz_reference", {"precursor","peptide"});

    defaults_.setValue("ignore_charge", "false", "For feature/consensus maps: Assign an ID independently of whether its charge state matches that of the (consensus) feature.");
    defaults_.setValidStrings("ignore_charge", {"true","false"});
    defaults_.setValue("match_meta_value", "", "For feature/consensus maps mapped by position (not TMT/iTRAQ data mapped by spectrum reference): "
                       "Assign an ID only to (consensus) features with the same value of this meta value, e.g. 'FAIMS_CV' to keep the "
                       "compensation voltages of FAIMS data apart. A missing value matches only a missing value. Empty: no restriction.", {"advanced"});

    defaultsToParam_();
  }

  IDMapper::IDMapper(const IDMapper& cp) :
    DefaultParamHandler(cp),
    rt_tolerance_(cp.rt_tolerance_),
    mz_tolerance_(cp.mz_tolerance_),
    measure_(cp.measure_),
    ignore_charge_(cp.ignore_charge_),
    match_meta_value_(cp.match_meta_value_)
  {
    updateMembers_();
  }

  IDMapper& IDMapper::operator=(const IDMapper& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    rt_tolerance_ = rhs.rt_tolerance_;
    mz_tolerance_ = rhs.mz_tolerance_;
    measure_ = rhs.measure_;
    ignore_charge_ = rhs.ignore_charge_;
    match_meta_value_ = rhs.match_meta_value_;
    updateMembers_();

    return *this;
  }

  void IDMapper::updateMembers_()
  {
    rt_tolerance_ = param_.getValue("rt_tolerance");
    mz_tolerance_ = param_.getValue("mz_tolerance");
    measure_ = param_.getValue("mz_measure") == "ppm" ? MEASURE_PPM : MEASURE_DA;
    ignore_charge_ = param_.getValue("ignore_charge") == "true";
    match_meta_value_ = param_.getValue("match_meta_value").toString();
  }

  void IDMapper::addIdentificationDataProcessing_(std::vector<DataProcessing>& data_processing, const std::vector<ProteinIdentification>& protein_ids)
  {
    for (const auto& prot_id : protein_ids)
    {
      DataProcessing dp;
      dp.getSoftware().setName(prot_id.getSearchEngine());
      dp.getSoftware().setVersion(prot_id.getSearchEngineVersion());
      dp.setCompletionTime(prot_id.getDateTime());
      dp.getProcessingActions().insert(DataProcessing::IDENTIFICATION);
      const auto& search_params = prot_id.getSearchParameters();
      if (!search_params.db.empty())
      {
        dp.setMetaValue("parameter: db", search_params.db);
      }
      if (!search_params.db_version.empty())
      {
        dp.setMetaValue("parameter: db_version", search_params.db_version);
      }
      data_processing.push_back(dp);
    }
  }

  void IDMapper::addIdentificationDataProcessing_(std::vector<DataProcessing>& data_processing, const IdentificationData& ids)
  {
    // one per legacy protein run, as for the protein identifications that export writes
    std::vector<ProteinIdentification> protein_ids;
    std::set<std::string> identifiers;
    for (const auto& run : ids.getRuns())
    {
      if (identifiers.insert(IdentificationDataAdapter::legacyIdentifier(run)).second)
      {
        protein_ids.push_back(IdentificationDataAdapter::settingsToLegacy(run));
      }
    }
    addIdentificationDataProcessing_(data_processing, protein_ids);
  }

  void IDMapper::annotate(AnnotatedMSRun& map,
    const PeptideIdentificationList& peptide_ids,
    const vector<ProteinIdentification>& protein_ids,
    const bool clear_ids,
    const bool map_ms1)
  {
    checkHits_(peptide_ids);
    SpectrumLookup lookup;

    if (clear_ids)
    { // start with empty IDs
      map.getPeptideIdentifications().clear();
      map.getProteinIdentifications().clear();
    }

    if (peptide_ids.empty()) return;

    // append protein identifications
    map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

    // AnnotatedMSRun will have one PeptideIdentification per spectrum (including ones without hits)
    map.getPeptideIdentifications().resize(map.getMSExperiment().getSpectra().size());
    
    // set up the lookup table for the spectra
    lookup.readSpectra(map.getMSExperiment());

    // remember which peptides were mapped (for stats later)
    unordered_set<Size> peptides_mapped;

    // store mapping of identification RT to index (ignore empty hits)
    multimap<double, Size> identifications_precursors;
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      if (peptide_ids[i].empty()) continue;      
      // mapping is done by either native id or by comparing peptide_id RT with experiment RT
      if (!peptide_ids[i].metaValueExists(Constants::UserParam::SPECTRUM_REFERENCE)) 
      { // use RT for mapping 
        identifications_precursors.insert(make_pair(peptide_ids[i].getRT(), i));
      } 
      else 
      { // use native id for mapping
        DataValue native_id = peptide_ids[i].getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE);
        try
        { // spectrum can be retrieved
          Size spectrum_idx = lookup.findByNativeID(native_id);
          // Since we now have only one PeptideIdentification per spectrum, we need to merge the hits
          PeptideIdentification& existing_id = map.getPeptideIdentifications()[spectrum_idx];
          existing_id.getHits().insert(existing_id.getHits().end(),
                                      peptide_ids[i].getHits().begin(),
                                      peptide_ids[i].getHits().end());
          peptides_mapped.insert(i);
        }
        catch (const Exception::ElementNotFound& /*e*/)
        { // use RT for mapping
          identifications_precursors.insert(make_pair(peptide_ids[i].getRT(), i));
        }
      }
    }

    if (!identifications_precursors.empty()) 
    {
      // store mapping of scan RT to index
      multimap<double, Size> experiment_precursors;
      for (Size i = 0; i < map.getMSExperiment().size(); i++)
      {
        experiment_precursors.insert(make_pair(map.getMSExperiment()[i].getRT(), i));
      }

      // note that mappings are sorted by key via multimap (we rely on that down below)

      // calculate the actual mapping
      multimap<double, Size>::const_iterator experiment_iterator = experiment_precursors.begin();
      multimap<double, Size>::const_iterator identifications_iterator = identifications_precursors.begin();
      // to achieve O(n) complexity we now move along the spectra
      // and for each spectrum we look at the peptide id's with the allowed RT range
      // once we finish a spectrum, we simply move back in the peptide id window a little to get from the
      // right end of the old interval to the left end of the new interval
      while (experiment_iterator != experiment_precursors.end())
      {
        // maybe we hit end() of IDs during the last scan - go back to a real value
        if (identifications_iterator == identifications_precursors.end())
        {
          --identifications_iterator; // this is valid, since we have at least one peptide ID
        }

        // go to left border of RT interval
        while (identifications_iterator != identifications_precursors.begin() &&
              (experiment_iterator->first - identifications_iterator->first) < rt_tolerance_) // do NOT use fabs() here, since we want the LEFT border
        {
          --identifications_iterator;
        }
        // ... we might have stepped too far left
        if (identifications_iterator != identifications_precursors.end() && ((experiment_iterator->first - identifications_iterator->first) > rt_tolerance_))
        {
          ++identifications_iterator; // get into interval again (we can potentially be at end() afterwards)
        }

        if (identifications_iterator == identifications_precursors.end())
        { // no more ID's, so we don't have any chance of matching the next spectra
          break; // ... do NOT put this block below, since hitting the end of ID's for one spec, still allows to match stuff in the next (when going to left border)
        }

        // run through RT interval
        while (identifications_iterator != identifications_precursors.end() &&
              (identifications_iterator->first - experiment_iterator->first) < rt_tolerance_) // fabs() not required here, since are definitely within left border, and wait until exceeding the right
        {
          bool success = map_ms1;
          if (!success)
          {
            for (const auto& precursor : map.getMSExperiment()[experiment_iterator->second].getPrecursors())
            {
              if (isMatch_(0, peptide_ids[identifications_iterator->second].getMZ(), precursor.getMZ()))
              {
                success = true;
                break;
              }
            }
          }

          if (success)
          {
            // Since we have only one PeptideIdentification per spectrum, we need to merge the hits
            PeptideIdentification& existing_id = map.getPeptideIdentifications()[experiment_iterator->second];
            existing_id.getHits().insert(existing_id.getHits().end(),
                                        peptide_ids[identifications_iterator->second].getHits().begin(),
                                        peptide_ids[identifications_iterator->second].getHits().end());
            peptides_mapped.insert(identifications_iterator->second);
          }
          ++identifications_iterator;
        }
        // we are at the right border now (or likely even beyond)
        ++experiment_iterator;
      }
    }

    // some statistics output
    OPENMS_LOG_INFO << "Peptides assigned to a precursor: " << peptides_mapped.size() << "\n"
             << "             Unassigned peptides: " << peptide_ids.size() - peptides_mapped.size() << "\n"
             << "       Unmapped (empty) peptides: " << peptide_ids.size() - identifications_precursors.size() << endl;
  }

  void IDMapper::annotate(AnnotatedMSRun& map, const FeatureMap& fmap, const bool clear_ids, const bool map_ms1)
  {
    // The identifications of a map with identification data, as peptide identifications
    std::optional<FeatureMap> exported;
    if (!fmap.getIdentificationData().empty())
    {
      exported.emplace(fmap);
      IdentificationDataConverter::exportFeatureIDs(*exported);
    }
    const FeatureMap& features = exported ? *exported : fmap;
    const vector<ProteinIdentification>& protein_ids = features.getProteinIdentifications();
    PeptideIdentificationList peptide_ids;

    for (FeatureMap::const_iterator it = features.begin(); it != features.end(); ++it)
    {
      const PeptideIdentificationList& pi = it->getPeptideIdentifications();
      for (PeptideIdentificationList::const_iterator itp = pi.begin(); itp != pi.end(); ++itp)
      {
        peptide_ids.push_back(*itp);
        // if pepID has no m/z or RT, use the values of the feature
        if (!itp->hasMZ()) peptide_ids.back().setMZ(it->getMZ());
        if (!itp->hasRT()) peptide_ids.back().setRT(it->getRT());
      }

    }
    annotate(map, peptide_ids, protein_ids, clear_ids, map_ms1);
  }

  enum class NATIVE_ID_TYPE
  {
    UNKNOWN, MS2IDMS3TMT, MS2IDTMT
  };

  NATIVE_ID_TYPE checkTMTType(const ConsensusMap& map)
  {
    for (auto & cf : map)
    {
      // check if the native id of an identifying spectrum is annotated
      if (cf.metaValueExists("id_scan_id")) // identifying MS2 spectrum in MS3 TMT
      {
        return NATIVE_ID_TYPE::MS2IDMS3TMT;
      }
      else if (cf.metaValueExists("scan_id")) // identifying MS2 spectrum in standard TMT
      {
        return NATIVE_ID_TYPE::MS2IDTMT;
      }
    }
    return NATIVE_ID_TYPE::UNKNOWN;
  }

  namespace
  {
    using ID = IdentificationData;

    /// An identification to map, with its run and the source that holds it
    struct Entry
    {
      const ID::Run* run = nullptr;
      const ID::Source* source = nullptr;
      const ID::Identification* query = nullptr;
    };

    std::set<std::string> uuidsOf(const ID& data)
    {
      std::set<std::string> uuids;
      for (const auto& run : data.getRuns())
      {
        uuids.insert(run.getUuid());
      }
      return uuids;
    }

    /// The identifications of the runs @p uuids of @p data, in the order of their IDs, then runs: the order of the
    /// peptide identifications they were imported from
    std::vector<Entry> entriesOf(const ID& data, const std::set<std::string>& uuids)
    {
      std::vector<std::tuple<UInt64, Size, Entry>> ordered;
      Size position = 0;
      for (const auto& run : data.getRuns())
      {
        if (uuids.contains(run.getUuid()))
        {
          for (const auto& source : run.getSources())
          {
            for (const auto& query : source.identifications)
            {
              ordered.emplace_back(query.getId().value, position, Entry {&run, &source, &query});
            }
          }
        }
        ++position;
      }
      std::stable_sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) {
        return std::make_pair(std::get<0>(a), std::get<1>(a)) < std::make_pair(std::get<0>(b), std::get<1>(b));
      });
      std::vector<Entry> entries;
      entries.reserve(ordered.size());
      for (const auto& item : ordered)
      {
        entries.push_back(std::get<2>(item));
      }
      return entries;
    }

    ID::QueryReference referenceOf(const Entry& entry)
    {
      return {entry.run->getUuid(), entry.query->getId()};
    }

    /// Link @p feature to the identification of @p entry with its matches
    void link(BaseFeature& feature, const Entry& entry)
    {
      feature.addIDQuery(referenceOf(entry));
      for (const auto& match : entry.query->getMatches())
      {
        feature.addIDMatch({entry.run->getUuid(), match.getId()});
      }
    }

    /// The file of the identification of @p entry, as for legacy peptide identifications (IdentifierMSRunMapper::getPrimaryMSRunPath())
    std::string sourcePath(const Entry& entry)
    {
      if (!entry.source->file.path.empty()) return entry.source->file.path;
      // a file that is not known: the first file of the run, or the deprecated meta value
      const StringList files = IdentificationDataAdapter::legacyFiles(*entry.run);
      if (!files.empty()) return files.front();
      return StringUtils::toStr(entry.query->getMetaValue(Constants::UserParam::BASE_NAME, ""));
    }

    /// An ID for new identifications of @p data, greater than those of all its identifications, so that export lists them after those
    UInt64 nextQueryId(const ID& data)
    {
      UInt64 next = 1;
      for (const auto& run : data.getRuns())
      {
        next = std::max(next, run.getNextQueryId());
      }
      return next;
    }

    /**
      @brief Identifications of precursors without identification

      They are in the first run of the map (a new run "UNKNOWN_SEARCH_RUN_IDENTIFIER" if there is none), as legacy
      peptide identifications without score type (IdentificationDataAdapter::unscoredRun()).
    */
    class Precursors
    {
    public:
      Precursors(ID& data, const PeakMap& spectra) :
        data_(data), spectra_(spectra)
      {
        if (data.getRuns().empty())
        {
          // a search run is mandatory, so we create one
          auto& run = data.addRun("UNKNOWN_SEARCH_RUN_IDENTIFIER");
          ID::RunSettings settings;
          settings.date = DateTime::now();
          run.setSettings(settings);
        }
        run_ = data.getRuns().front().getIdentifier();
      }

      /// A new identification of precursor @p mz of spectrum @p spectrum_index
      ID::QueryReference add(Size spectrum_index, double mz, std::optional<UInt64> map_index = std::nullopt)
      {
        if (!next_)
        {
          run_ = IdentificationDataAdapter::unscoredRun(data_, run_).getIdentifier();
          next_ = nextQueryId(data_);
        }
        auto& run = data_.getRun(run_);
        ID::Observation observation;
        observation.rt = spectra_[spectrum_index].getRT();
        observation.mz = mz;
        observation.setMetaValue("spectrum_index", spectrum_index);
        observation.data_id = spectra_[spectrum_index].getNativeID();
        if (map_index)
        {
          // we use no underscore here to be compatible with linkers
          observation.setMetaValue("map_index", *map_index);
        }
        // the file of the precursor is not known (if the run has several), as for legacy peptide identifications without 'id_merge_index'
        const auto source = IdentificationDataAdapter::legacySource(run, IdentificationDataAdapter::legacyFiles(run).size(), PeptideIdentification());
        return {run.getUuid(), run.importIdentification(source, ID::QueryId {next_++}, std::move(observation))};
      }

    private:
      ID& data_;
      const PeakMap& spectra_;
      std::string run_;
      UInt64 next_ = 0;
    };

    /**
      @brief Annotate @p map as identification data with @p annotate

      A map with identification data is annotated as it is. Otherwise, its peptide identifications are moved into its
      identification data and back afterwards (also if @p annotate throws), so it holds peptide identifications.
    */
    template<class Map>
    void annotateAsIdentificationData(Map& map, const std::function<void(Map&)>& annotate)
    {
      if (!map.getIdentificationData().empty() && !IdentificationDataConverter::hasPeptideIdentifications(map))
      {
        annotate(map);
        return;
      }
      IdentificationDataConverter::moveToIdentificationData(map);
      const auto restore = [&map]() {
        if constexpr (std::is_same_v<Map, FeatureMap>) IdentificationDataConverter::exportFeatureIDs(map);
        else IdentificationDataConverter::exportConsensusIDs(map);
      };
      try
      {
        annotate(map);
      }
      catch (...)
      {
        restore();
        throw;
      }
      restore();
    }

    template<class Map>
    void checkNative(const Map& map)
    {
      if (IdentificationDataConverter::hasPeptideIdentifications(map))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "IDMapper: the map has peptide identifications; annotating it with identification data needs them as "
                                          "identification data (see IdentificationDataConverter::moveToIdentificationData())");
      }
    }
  } // namespace

  void IDMapper::annotate(
    ConsensusMap& map,
    const PeptideIdentificationList& ids,
    const vector<ProteinIdentification>& protein_ids,
    bool measure_from_subelements,
    bool annotate_ids_with_subelements,
    const PeakMap& spectra)
  {
    // validate "RT" and "MZ" metavalues exist
    checkHits_(ids);
    auto imported = IdentificationDataAdapter::fromLegacy(protein_ids, ids);

    // preserve data processing from identification runs (search engine, database, etc.)
    addIdentificationDataProcessing_(map.getDataProcessing(), protein_ids);

    annotateAsIdentificationData<ConsensusMap>(map, [&](ConsensusMap& native) {
      annotate_(native, std::move(imported), measure_from_subelements, annotate_ids_with_subelements, spectra);
    });
  }

  void IDMapper::annotate(ConsensusMap& map, const IdentificationData& ids, bool measure_from_subelements, bool annotate_ids_with_subelements, const PeakMap& spectra)
  {
    checkNative(map);
    annotate_(map, ids, measure_from_subelements, annotate_ids_with_subelements, spectra);

    // preserve data processing from identification runs (search engine, database, etc.)
    addIdentificationDataProcessing_(map.getDataProcessing(), ids);
  }

  void IDMapper::annotate_(ConsensusMap& map, IdentificationData ids, bool measure_from_subelements, bool annotate_ids_with_subelements, const PeakMap& spectra)
  {
    checkHits_(ids);
    const vector<Size> unidentified = mapPrecursorsToIdentifications(spectra, ids).unidentified;
    // whether spectrum references of TMT/iTRAQ data are scan numbers is told by the first identification
    std::string first_reference;
    if (const auto all = entriesOf(ids, uuidsOf(ids)); !all.empty())
    {
      first_reference = all.front().query->data_id;
    }

    // identifications without matches are not mapped
    ids.eraseIdentifications([](const ID::Run&, const ID::Identification& query) { return query.getMatches().empty(); });
    const std::set<std::string> uuids = uuidsOf(ids);
    ID& data = map.getIdentificationData();
    data.merge(ids);
    const std::vector<Entry> entries = entriesOf(data, uuids);

    // keep track of assigned/unassigned identifications.
    // maps entry index to number of assignments to a feature
    std::unordered_map<Size, Size> assigned_ids;

    // keep track of assigned/unassigned precursors
    std::unordered_map<Size, Size> assigned_precursors;

    // store which identifications fit which feature (and avoid double entries)
    // consensusMap -> {entry index}
    vector<set<size_t>> mapping(map.size());

    DoubleList mz_values;
    double rt_pep;
    IntList charges;

    // for statistics
    Size id_matches_none(0), id_matches_single(0), id_matches_multiple(0);

    NATIVE_ID_TYPE native_id_type = checkTMTType(map);

    // We have TMT data: spectrum references annotated at consensus feature and in id
    // We can directly map by native id
    if ((native_id_type != NATIVE_ID_TYPE::UNKNOWN) )
    {
      // build map from file to identification
      std::map<std::string, std::unordered_map<std::string, const Entry*>> file2nativeid2entry;
      bool has_spectrum_references{false};
      bool lookForScanNrsAsIntegers = false;

      for (const Entry& entry : entries)
      {
        const std::string& spectrum_reference = entry.query->data_id;
        // missing file origin is fine, but we need a spectrum_reference if we want to build the map
        if (spectrum_reference.empty()) continue;
        // TODO make a unique decision in the whole class on if to extract by scan number or full string?
        if (!lookForScanNrsAsIntegers)
        {
          // check if spectrum reference is a string that just contains a number
          try
          {
            StringUtils::toInt64(first_reference);
            lookForScanNrsAsIntegers = true;
          }
          catch (...)
          {
            lookForScanNrsAsIntegers = false;
          }
        }
        auto& inner_map = file2nativeid2entry[File::basename(sourcePath(entry))];
        auto result = inner_map.insert({spectrum_reference, &entry});
        if (!result.second)
        {
          OPENMS_LOG_WARN << "Duplicate spectrum reference detected: "<< spectrum_reference << "\n";
        }
        has_spectrum_references = true;
      }

      if (!has_spectrum_references)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No spectrum references in ID file used in id mapping TMT/iTRAQ data.");
      }

      if (measure_from_subelements)
      {
        OPENMS_LOG_WARN << "IDMapper is configured to measure from subelements. Because the data looks like TMT/iTRAQ this option will be ignored." << std::endl;
      }

      if (!ignore_charge_)
      {
        OPENMS_LOG_WARN << "IDMapper is configured to validate charges. Because the data looks like TMT/iTRAQ this option will be ignored."  << std::endl;
      }

      // the identifications that a consensus feature links
      std::set<ID::QueryReference> linked;
      // Default-constructed, so empty() holds until the first scan_id sets it below. A compiled
      // empty pattern ("") is not empty(): with it the fallback never ran, and extractScanNumber()
      // found no capture group and threw. Declared outside the loop, the regex is derived once.
      RegularExpression scanregex;
      for (auto& cf : map)
      {
        const auto first_channel = *cf.getFeatures().begin();
        std::string filename = File::basename(map.getColumnHeaders()[first_channel.getMapIndex()].filename); // all channels are associated with same file in TMT/iTRAQ

        std::string cf_scan_id_key_name = (native_id_type == NATIVE_ID_TYPE::MS2IDMS3TMT) ? "id_scan_id" : "scan_id";
        std::string cf_scan_id = StringUtils::toStr(cf.getMetaValue(cf_scan_id_key_name, ""));
        if (!cf_scan_id.empty())
        {
          // This assumes all scan_ids are of the same structure
          if (lookForScanNrsAsIntegers && scanregex.empty()) { scanregex.assign(SpectrumLookup::getRegExFromNativeID(cf_scan_id)); }
          if (auto run_it = file2nativeid2entry.find(filename); run_it != file2nativeid2entry.end()) // TMT/iTRAQ run has identifications
          {
            if (auto scanid_it = run_it->second.find(cf_scan_id); scanid_it != run_it->second.end()) // TMT/iTRAQ run has scan_id with identification
            {
              link(cf, *scanid_it->second);
              linked.insert(referenceOf(*scanid_it->second));
              ++id_matches_single; // in TMT we only match to single consensus feature
            }
            // look for only the scan_number in case the search engine only extracted this (e.g. Sage)
            else if (lookForScanNrsAsIntegers)
            {
              // A WIFF native ID ("sample=1 period=1 cycle=96 experiment=1") holds two numbers, and
              // the generic regex would take the last, the experiment. Its scan number is
              // cycle * 1000 + experiment (96001), which the accession-based overload computes.
              const Int scan_number = StringUtils::hasSubstring(cf_scan_id, "cycle=")
                ? SpectrumLookup::extractScanNumber(cf_scan_id, "MS:1000770")
                : SpectrumLookup::extractScanNumber(cf_scan_id, scanregex, false);
              auto scanid_it = run_it->second.find(StringUtils::toStr(scan_number));
              if(scanid_it != run_it->second.end())
              {
                link(cf, *scanid_it->second);
                linked.insert(referenceOf(*scanid_it->second));
                ++id_matches_single; // in TMT we only match to single consensus feature
              }
            }
          } // else identification file does not contained scan id (e.g. was removed)
          else
          {
            OPENMS_LOG_WARN << "ConsensusMap for TMT/iTRAQ experiment contains scan identifier '" << cf_scan_id
                          << "' quantified in file '" << filename
                          << "' but there is no matching identification."
                          << std::endl;
          }
        }
        else // missing spectrum id annotation
        {
            OPENMS_LOG_WARN << "ConsensusMap for TMT/iTRAQ experiment is missing the scan identifier meta value '" << cf_scan_id << "'"
                          << std::endl;
        }
      }
      // TMT/iTRAQ data has no unassigned identifications: those that map to no consensus feature are not kept
      data.eraseIdentifications([&](const ID::Run& run, const ID::Identification& query) {
        return uuids.contains(run.getUuid()) && !linked.contains({run.getUuid(), query.getId()});
      });
    }
    else
    { // non TMT data (e.g., label-free)
      // assignments that are identifications of their own, with the map index of the matching subelement
      struct Assignment
      {
        Size feature;
        Size entry;
        UInt64 map_index;
      };
      std::vector<Assignment> own_assignments;

      for (Size i = 0; i < entries.size(); ++i)
      {
        getIDDetails_(*entries[i].run, *entries[i].query, rt_pep, mz_values, charges);

        bool id_mapped(false);

        // iterate over the features
        for (Size cm_index = 0; cm_index < map.size(); ++cm_index)
        {
          if (!sameMetaValue_(map[cm_index], *entries[i].query)) continue;

          // if set to TRUE, we leave the i_mz-loop as we added the whole ID with all hits
          bool was_added = false; // was current pep-m/z matched?!

          // iterate over m/z values of pepIds
          for (Size i_mz = 0; i_mz < mz_values.size(); ++i_mz)
          {
            double mz_pep = mz_values[i_mz];

            // charge states to use for checking:
            IntList current_charges;
            if (!ignore_charge_)
            {
              // if "mz_ref." is "precursor", we have only one m/z value to check,
              // but still one charge state per peptide hit that could match:
              if (mz_values.size() == 1)
              {
                current_charges = charges;
              }
              else
              {
                current_charges.push_back(charges[i_mz]);
              }
              current_charges.push_back(0); // "not specified" always matches
            }

            //check if we compare distance from centroid or subelements
            if (!measure_from_subelements)
            {
              if (
                  isMatch_(rt_pep - map[cm_index].getRT(), mz_pep, map[cm_index].getMZ()) &&
                  (ignore_charge_ || ListUtils::contains(current_charges, map[cm_index].getCharge()))
                  )
              {
                id_mapped = true;
                was_added = true;
                link(map[cm_index], entries[i]);
                ++assigned_ids[i];
              }
            }
            else
            {
              for (ConsensusFeature::HandleSetType::const_iterator it_handle = map[cm_index].getFeatures().begin();
                  it_handle != map[cm_index].getFeatures().end();
                  ++it_handle)
              {
                if (isMatch_(rt_pep - it_handle->getRT(), mz_pep, it_handle->getMZ()) &&
                    (ignore_charge_ || ListUtils::contains(current_charges, it_handle->getCharge())))
                {
                  id_mapped = true;
                  was_added = true;
                  if (!mapping[cm_index].contains(i))
                  {
                    if (annotate_ids_with_subelements)
                    {
                      // the identification gets the map index of the matching subelement, so this assignment is an identification of its own
                      own_assignments.push_back({cm_index, i, it_handle->getMapIndex()});
                    }
                    else
                    {
                      link(map[cm_index], entries[i]);
                    }
                    ++assigned_ids[i];
                    mapping[cm_index].insert(i);
                  }
                  break; // we added this peptide already.. no need to check other handles
                }
              }
              // continue to here
            }

            if (was_added) break;

          } // m/z values to check

          // break to here

        } // features

        // the id has not been mapped to any consensus feature, so it stays unassigned
        if (!id_mapped)
        {
          ++id_matches_none;
        }
      } // Identifications

      for (auto aid : assigned_ids)
      {
        if (aid.second == 1)
        {
          ++id_matches_single;
        }
        else if (aid.second > 1)
        {
          ++id_matches_multiple;
        }
      }

      if (!own_assignments.empty())
      {
        // Copies of the assigned identifications (with their matches and scores) replace them, in the order of assignment.
        struct Copy
        {
          std::string run;
          ID::SourceId source;
          ID::Observation observation;
          std::vector<std::pair<ID::MatchData, std::vector<std::optional<double>>>> matches;
          Size feature;
        };
        std::vector<Copy> copies;
        std::set<ID::QueryReference> copied;
        for (const auto& assignment : own_assignments)
        {
          const Entry& entry = entries[assignment.entry];
          Copy copy {entry.run->getIdentifier(), entry.source->id, entry.query->getObservation(), {}, assignment.feature};
          // we use no underscore here to be compatible with linkers
          copy.observation.setMetaValue("map_index", assignment.map_index);
          for (const auto& match : entry.query->getMatches())
          {
            copy.matches.emplace_back(match.getData(), entry.run->getScores(match));
          }
          copies.push_back(std::move(copy));
          copied.insert(referenceOf(entry));
        }
        UInt64 next = nextQueryId(data);
        for (auto& copy : copies)
        {
          auto& run = data.getRun(copy.run);
          const auto query = run.importIdentification(copy.source, ID::QueryId {next++}, std::move(copy.observation));
          auto& feature = map[copy.feature];
          feature.addIDQuery({run.getUuid(), query});
          for (const auto& [match, scores] : copy.matches)
          {
            feature.addIDMatch({run.getUuid(), run.addMatch(query, match, scores)});
          }
        }
        data.eraseIdentifications([&](const ID::Run& run, const ID::Identification& query) { return copied.contains({run.getUuid(), query.getId()}); });
      }
    }

    if (!entries.empty() && !spectra.empty())
    {
      OPENMS_LOG_INFO << "Mapping " << entries.size() << "PeptideIdentifications to " << spectra.size() << " spectra." << endl;

      const auto state = mapPrecursorsToIdentifications(spectra, ids);
      OPENMS_LOG_INFO << "Identification state of spectra: \n"
               << "Unidentified: " << state.unidentified.size() << "\n"
               << "Identified:   " << state.identified.size() << "\n"
               << "No precursor: " << state.no_precursors.size() << endl;
    }

    // identifications of unidentified precursors need a search run
    std::optional<Precursors> precursors;
    if (!unidentified.empty())
    {
      precursors.emplace(data, spectra);
    }

    // for statistics:
    Size spectrum_matches_none(0), spectrum_matches_single(0), spectrum_matches_multiple(0);

    // are there any mapped but unidentified precursors?
    for (Size ui = 0; ui != unidentified.size(); ++ui)
    {
      Size spectrum_index = unidentified[ui];
      const MSSpectrum& spectrum = spectra[spectrum_index];
      const vector<Precursor>& precursor_list = spectrum.getPrecursors();

      bool precursor_mapped(false);

      // check if precursor has been identified
      for (Size i_p = 0; i_p < precursor_list.size(); ++i_p)
      {
        // check by precursor mass and spectrum RT
        double mz_p = precursor_list[i_p].getMZ();
        int z_p = precursor_list[i_p].getCharge();
        double rt_value = spectrum.getRT();

        // the identification of the precursor that the matching consensus features link
        std::optional<ID::QueryReference> precursor_id;

        // iterate over the consensus features
        for (Size cm_index = 0; cm_index < map.size(); ++cm_index)
        {
          // charge states to use for checking:
          IntList current_charges;
          if (!ignore_charge_)
          {
            current_charges.push_back(z_p);
            current_charges.push_back(0); // "not specified" always matches
          }

          // check if we compare distance from centroid or subelements
          if (!measure_from_subelements) // measure from centroid
          {
            if (isMatch_(rt_value - map[cm_index].getRT(), mz_p, map[cm_index].getMZ()) && (ignore_charge_ || ListUtils::contains(current_charges, map[cm_index].getCharge())))
            {
              if (!precursor_id) precursor_id = precursors->add(spectrum_index, mz_p);
              map[cm_index].addIDQuery(*precursor_id);
              ++assigned_precursors[spectrum_index];
              precursor_mapped = true;
            }
          }
          else // measure from subelements
          {
            for (ConsensusFeature::HandleSetType::const_iterator it_handle = map[cm_index].getFeatures().begin();
                 it_handle != map[cm_index].getFeatures().end();
                 ++it_handle)
            {
              if (isMatch_(rt_value - it_handle->getRT(), mz_p, it_handle->getMZ())  && (ignore_charge_ || ListUtils::contains(current_charges, it_handle->getCharge())))
              {
                if (annotate_ids_with_subelements)
                {
                  // store the map index the precursor was mapped to: an identification for this subelement
                  map[cm_index].addIDQuery(precursors->add(spectrum_index, mz_p, it_handle->getMapIndex()));
                }
                else
                {
                  if (!precursor_id) precursor_id = precursors->add(spectrum_index, mz_p);
                  map[cm_index].addIDQuery(*precursor_id);
                }
                ++assigned_precursors[spectrum_index];
                precursor_mapped = true;
              }
            }
          }
        } // m/z values to check
      }
      if (!precursor_mapped) ++spectrum_matches_none;
    }

    for (auto apc : assigned_precursors)
    {
      if (apc.second == 1)
      {
        ++spectrum_matches_single;
      }
      else if (apc.second > 1)
      {
        ++spectrum_matches_multiple;
      }
    }

    // some statistics output
    if (!entries.empty())
    {
      OPENMS_LOG_INFO << "Unassigned peptides: " << id_matches_none << "\n"
               << "Peptides assigned to exactly one feature: " << id_matches_single << "\n"
               << "Peptides assigned to multiple features: " << id_matches_multiple << "\n";
    }

    if (!spectra.empty())
    {
      OPENMS_LOG_INFO << "Unassigned precursors without identification: " << spectrum_matches_none << "\n"
               << "Unidentified precursor assigned to exactly one feature: " << spectrum_matches_single << "\n"
               << "Unidentified precursor assigned to multiple features: " << spectrum_matches_multiple << endl;
    }
  }

  void IDMapper::annotate(FeatureMap& map,
    const PeptideIdentificationList& ids,
    const vector<ProteinIdentification>& protein_ids,
    bool use_centroid_rt,
    bool use_centroid_mz,
    const PeakMap& spectra)
  {
    checkHits_(ids); // check RT and m/z are present
    auto imported = IdentificationDataAdapter::fromLegacy(protein_ids, ids);

    // preserve data processing from identification runs (search engine, database, etc.)
    addIdentificationDataProcessing_(map.getDataProcessing(), protein_ids);

    annotateAsIdentificationData<FeatureMap>(map, [&](FeatureMap& native) {
      annotate_(native, std::move(imported), use_centroid_rt, use_centroid_mz, spectra);
    });
  }

  void IDMapper::annotate(FeatureMap& map, const IdentificationData& ids, bool use_centroid_rt, bool use_centroid_mz, const PeakMap& spectra)
  {
    checkNative(map);
    annotate_(map, ids, use_centroid_rt, use_centroid_mz, spectra);

    // preserve data processing from identification runs (search engine, database, etc.)
    addIdentificationDataProcessing_(map.getDataProcessing(), ids);
  }

  void IDMapper::annotate_(FeatureMap& map, IdentificationData ids, bool use_centroid_rt, bool use_centroid_mz, const PeakMap& spectra)
  {
    checkHits_(ids); // check RT and m/z are present
    const vector<Size> unidentified = mapPrecursorsToIdentifications(spectra, ids).unidentified;

    // identifications without matches are not mapped
    ids.eraseIdentifications([](const ID::Run&, const ID::Identification& query) { return query.getMatches().empty(); });
    const std::set<std::string> uuids = uuidsOf(ids);
    ID& data = map.getIdentificationData();
    data.merge(ids);
    const std::vector<Entry> entries = entriesOf(data, uuids);

    // check if all features have at least one convex hull
    // if not, use the centroid and the given tolerances
    if (!(use_centroid_rt && use_centroid_mz))
    {
      for (Feature& f_it : map)
      {
        if (f_it.getConvexHulls().empty())
        {
          use_centroid_rt = true;
          use_centroid_mz = true;
          OPENMS_LOG_WARN << "IDMapper warning: at least one feature has no convex hull - using centroid coordinates for matching" << endl;
          break;
        }
      }
    }

    bool use_avg_mass = false;           // use avg. peptide masses for matching?
    if (use_centroid_mz && (param_.getValue("mz_reference") == "peptide"))
    {
      // if possible, check which m/z value is reported for features,
      // so the appropriate peptide mass can be used for matching
      use_avg_mass = checkMassType_(map.getDataProcessing());
    }

    // calculate feature bounding boxes only once:
    vector<DBoundingBox<2> > boxes;
    double min_rt = numeric_limits<double>::max();
    double max_rt = -numeric_limits<double>::max();
    // cout << "Precomputing bounding boxes..." << endl;
    boxes.reserve(map.size());
    for (Feature& f_it : map)
    {
      DBoundingBox<2> box;
      if (!(use_centroid_rt && use_centroid_mz))
      {
        box = f_it.getConvexHull().getBoundingBox();
      }
      if (use_centroid_rt)
      {
        box.setMinX(f_it.getRT());
        box.setMaxX(f_it.getRT());
      }
      if (use_centroid_mz)
      {
        box.setMinY(f_it.getMZ());
        box.setMaxY(f_it.getMZ());
      }
      increaseBoundingBox_(box);
      boxes.push_back(box);

      min_rt = min(min_rt, box.minPosition().getX());
      max_rt = max(max_rt, box.maxPosition().getX());
    }

    // hash bounding boxes of features by RT:
    // RT range is partitioned into slices (bins) of 1 second; every feature
    // that overlaps a certain slice is hashed into the corresponding bin
    vector<vector<SignedSize> > hash_table;
    // make sure the RT hash table has indices >= 0 and doesn't waste space
    // in the beginning:
    SignedSize offset(0);

    if (!map.empty())
    {
      // cout << "Setting up hash table..." << endl;
      offset = SignedSize(floor(min_rt));
      // this only works if features were found
      hash_table.resize(SignedSize(floor(max_rt)) - offset + 1);
      for (Size index = 0; index < boxes.size(); ++index)
      {
        const DBoundingBox<2> & box = boxes[index];
        for (SignedSize i = SignedSize(floor(box.minPosition().getX()));
             i <= SignedSize(floor(box.maxPosition().getX())); ++i)
        {
          hash_table[i - offset].push_back(index);
        }
      }
    }
    else
    {
      OPENMS_LOG_WARN << "IDMapper received an empty FeatureMap! All peptides are mapped as 'unassigned'!" << endl;
    }

    // for statistics:
    Size matches_none = 0, matches_single = 0, matches_multi = 0;

    // cout << "Finding matches..." << endl;
    // iterate over identifications (an identification that matches no feature stays unassigned):
    for (const Entry& entry : entries)
    {
      DoubleList mz_values;
      double rt_value;
      IntList charges;
      getIDDetails_(*entry.run, *entry.query, rt_value, mz_values, charges, use_avg_mass);

      if ((rt_value < min_rt) || (rt_value > max_rt)) // RT out of bounds
      {
        ++matches_none;
        continue;
      }

      // iterate over candidate features:
      Size index = SignedSize(floor(rt_value)) - offset;
      Size matching_features = 0;
      for (SignedSize& hash_it : hash_table[index])
      {
        Feature & feat = map[hash_it];
        if (!sameMetaValue_(feat, *entry.query)) continue;

        // need to check the charge state?
        bool check_charge = !ignore_charge_;
        if (check_charge && (mz_values.size() == 1))               // check now
        {
          if (!ListUtils::contains(charges, feat.getCharge())) continue;
          check_charge = false;                 // don't need to check later
        }

        // iterate over m/z values (only one if "mz_ref." is "precursor"):
        Size l_index = 0;
        for (DoubleList::iterator mz_it = mz_values.begin();
             mz_it != mz_values.end(); ++mz_it, ++l_index)
        {
          if (check_charge && (charges[l_index] != feat.getCharge()))
          {
            continue;                   // charge states need to match
          }

          DPosition<2> id_pos(rt_value, *mz_it);
          if (boxes[hash_it].encloses(id_pos))                 // potential match
          {
            if (use_centroid_mz)
            {
              // only one m/z value to check, which was already incorporated
              // into the overall bounding box -> success!
              link(feat, entry);
              ++matching_features;
              break;                     // "mz_it" loop
            }
            // else: check all the mass traces
            bool found_match = false;
            for (vector<ConvexHull2D>::iterator ch_it =
                 feat.getConvexHulls().begin(); ch_it !=
                 feat.getConvexHulls().end(); ++ch_it)
            {
              DBoundingBox<2> box = ch_it->getBoundingBox();
              if (use_centroid_rt)
              {
                box.setMinX(feat.getRT());
                box.setMaxX(feat.getRT());
              }
              increaseBoundingBox_(box);
              if (box.encloses(id_pos)) // success!
              {
                link(feat, entry);
                ++matching_features;
                found_match = true;
                break; // "ch_it" loop
              }
            }
            if (found_match) break; // "mz_it" loop
          }
        }
      }
      if (matching_features == 0)
      {
        ++matches_none;
      }
      else if (matching_features == 1)
      {
        ++matches_single;
      }
      else
      {
        ++matches_multi;
      }
    }

    // map all unidentified precursor to features
    Size spectrum_matches_none(0);
    Size spectrum_matches_single(0);
    Size spectrum_matches_multi(0);

    // identifications of unidentified precursors need a search run
    std::optional<Precursors> precursors;
    if (!unidentified.empty())
    {
      precursors.emplace(data, spectra);
    }

    // are there any mapped but unidentified precursors?
    for (Size i = 0; i != unidentified.size(); ++i)
    {
      Size spectrum_index = unidentified[i];
      const MSSpectrum& spectrum = spectra[spectrum_index];
      const vector<Precursor>& precursor_list = spectrum.getPrecursors();

      // check if precursor has been identified
      for (Size i_p = 0; i_p < precursor_list.size(); ++i_p)
      {
        // check by precursor mass and spectrum RT
        double mz_p = precursor_list[i_p].getMZ();
        double rt_value = spectrum.getRT();
        int z_p = precursor_list[i_p].getCharge();

        if ((rt_value < min_rt) || (rt_value > max_rt)) // RT out of bounds
        {
          ++spectrum_matches_none;
          continue;
        }

        // iterate over candidate features:
        Size index = SignedSize(floor(rt_value)) - offset;
        Size matching_features = 0;

        for (SignedSize& hash_it : hash_table[index])
        {
          Feature & feat = map[hash_it];

          // (optionally) check charge state
          if (!ignore_charge_)
          {
            if (std::abs(z_p) != std::abs(feat.getCharge())) continue;
          }

          DPosition<2> id_pos(rt_value, mz_p);

          if (boxes[hash_it].encloses(id_pos)) // potential match
          {
            if (use_centroid_mz)
            {
              // only one m/z value to check, which was already incorporated
              // into the overall bounding box -> success!
              feat.addIDQuery(precursors->add(spectrum_index, mz_p));
              ++matching_features;
              break; // "mz_it" loop
            }
            // else: check all the mass traces
            bool found_match = false;
            for (vector<ConvexHull2D>::iterator ch_it =
                  feat.getConvexHulls().begin(); ch_it !=
                  feat.getConvexHulls().end(); ++ch_it)
            {
              DBoundingBox<2> box = ch_it->getBoundingBox();
              if (use_centroid_rt)
              {
                box.setMinX(feat.getRT());
                box.setMaxX(feat.getRT());
              }
              increaseBoundingBox_(box);
              if (box.encloses(id_pos)) // success!
              {
                feat.addIDQuery(precursors->add(spectrum_index, mz_p));
                ++matching_features;
                found_match = true;
                break; // "ch_it" loop
              }
            }

            if (found_match) break; // "mz_it" loop
          }
        }

        if (matching_features == 0)
        {
          ++spectrum_matches_none;
        }
        else if (matching_features == 1)
        {
          ++spectrum_matches_single;
        }
        else
        {
          ++spectrum_matches_multi;
        }
      }
    }

    // some statistics output
    OPENMS_LOG_INFO << "Unassigned peptides: " << matches_none << "\n"
    << "Peptides assigned to exactly one feature: " << matches_single << "\n"
    << "Peptides assigned to multiple features: " << matches_multi << "\n";

    OPENMS_LOG_INFO << "Unassigned and unidentified precursors: " << spectrum_matches_none << "\n"
    << "Unidentified precursor assigned to exactly one feature: " << spectrum_matches_single << "\n"
    << "Unidentified precursor assigned to multiple features: " << spectrum_matches_multi << "\n";

    OPENMS_LOG_INFO << map.getAnnotationStatistics() << endl;
  }

  IDMapper::PeptideIdentificationListState IDMapper::mapPrecursorsToIdentifications(const PeakMap& spectra, const IdentificationData& ids, double mz_tol,
                                                                                    double rt_tol)
  {
    // positions (m/z, RT) of the identifications with matches: those without do not identify a spectrum
    std::vector<std::pair<double, double>> positions;
    for (const auto& run : ids.getRuns())
    {
      for (const auto& source : run.getSources())
      {
        for (const auto& query : source.identifications)
        {
          if (query.getMatches().empty()) continue;
          positions.emplace_back(query.mz.value_or(std::numeric_limits<double>::quiet_NaN()), query.rt.value_or(std::numeric_limits<double>::quiet_NaN()));
        }
      }
    }
    PeptideIdentificationListState ret;
    for (Size spectrum_index = 0; spectrum_index < spectra.size(); ++spectrum_index)
    {
      const MSSpectrum& spectrum = spectra[spectrum_index];
      if (spectrum.getPrecursors().empty())
      {
        ret.no_precursors.push_back(spectrum_index);
        continue;
      }
      // check if a precursor has been identified, by precursor mass and spectrum RT
      const bool identified = std::any_of(spectrum.getPrecursors().begin(), spectrum.getPrecursors().end(), [&](const Precursor& precursor) {
        return std::any_of(positions.begin(), positions.end(), [&](const auto& position) {
          return fabs(position.first - precursor.getMZ()) < mz_tol && fabs(spectrum.getRT() - position.second) < rt_tol;
        });
      });
      (identified ? ret.identified : ret.unidentified).push_back(spectrum_index);
    }
    return ret;
  }

  double IDMapper::getAbsoluteMZTolerance_(const double mz) const
  {
    if (measure_ == MEASURE_PPM)
    {
      return mz * mz_tolerance_ / 1e6;
    }
    else if (measure_ == MEASURE_DA)
    {
      return mz_tolerance_;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IDMapper::getAbsoluteTolerance_(): illegal internal state of measure_!",StringUtils::toStr(measure_));
  }

  bool IDMapper::isMatch_(const double rt_distance, const double mz_theoretical, const double mz_observed) const
  {
    if (measure_ == MEASURE_PPM)
    {
      return (fabs(rt_distance) <= rt_tolerance_) && (Math::getPPMAbs(mz_observed, mz_theoretical) <= mz_tolerance_);
    }
    else if (measure_ == MEASURE_DA)
    {
      return (fabs(rt_distance) <= rt_tolerance_) && (fabs(mz_theoretical - mz_observed) <= mz_tolerance_);
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IDMapper::getAbsoluteTolerance_(): illegal internal state of measure_!",StringUtils::toStr(measure_));
  }

  bool IDMapper::sameMetaValue_(const MetaInfoInterface& feature, const MetaInfoInterface& identification) const
  {
    if (match_meta_value_.empty()) return true;
    const bool feature_has = feature.metaValueExists(match_meta_value_);
    if (feature_has != identification.metaValueExists(match_meta_value_)) return false;
    return ! feature_has || feature.getMetaValue(match_meta_value_) == identification.getMetaValue(match_meta_value_);
  }

  void IDMapper::checkHits_(const PeptideIdentificationList& ids) const
  {
    for (Size i = 0; i < ids.size(); ++i)
    {
      if (!ids[i].hasRT())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IDMapper: 'RT' information missing for peptide identification!");
      }
      if (!ids[i].hasMZ())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IDMapper: 'MZ' information missing for peptide identification!");
      }
    }
  }

  void IDMapper::checkHits_(const IdentificationData& ids) const
  {
    for (const auto& run : ids.getRuns())
    {
      for (const auto& source : run.getSources())
      {
        for (const auto& query : source.identifications)
        {
          if (!query.rt)
          {
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IDMapper: 'RT' information missing for peptide identification!");
          }
          if (!query.mz)
          {
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IDMapper: 'MZ' information missing for peptide identification!");
          }
        }
      }
    }
  }

  void IDMapper::getIDDetails_(const IdentificationData::Run& run, const IdentificationData::Identification& id, double& rt_pep, DoubleList& mz_values,
                               IntList& charges, bool use_avg_mass) const
  {
    mz_values.clear();
    charges.clear();

    rt_pep = *id.rt;

    // collect m/z values of the identification
    if (param_.getValue("mz_reference") == "precursor") // use precursor m/z of the identification
    {
      mz_values.push_back(*id.mz);
    }

    for (const auto& match : id.getMatches())
    {
      Int charge = match.charge;
      charges.push_back(charge);

      if (param_.getValue("mz_reference") == "peptide") // use mass of each match (assuming H+ adducts)
      {
        const AASequence sequence = IdentificationDataAdapter::materializePeptide(run, match, *run.getPrimaryScore()).getSequence();
        double mass = use_avg_mass ?
                      sequence.getAverageWeight(Residue::Full, charge) :
                      sequence.getMonoWeight(Residue::Full, charge);

        mz_values.push_back(mass / (double) charge);
      }
    }
  }

  void IDMapper::increaseBoundingBox_(DBoundingBox<2>& box)
  {
    DPosition<2> sub_min(rt_tolerance_,
                         getAbsoluteMZTolerance_(box.minPosition().getY())),
    add_max(rt_tolerance_, getAbsoluteMZTolerance_(box.maxPosition().getY()));

    box.setMin(box.minPosition() - sub_min);
    box.setMax(box.maxPosition() + add_max);
  }

  bool IDMapper::checkMassType_(const vector<DataProcessing>& processing) const
  {
    bool use_avg_mass = false;
    std::string before;
    for (const DataProcessing& proc_it : processing)
    {
      if (proc_it.getSoftware().getName() == "FeatureFinder")
      {
        std::string reported_mz = StringUtils::toStr(proc_it.getMetaValue("parameter: algorithm:feature:reported_mz"));
        if (reported_mz.empty())
          continue; // parameter info not available
        if (!before.empty() && (reported_mz != before))
        {
          OPENMS_LOG_WARN << "The m/z values reported for features in the input seem to be of different types (e.g. monoisotopic/average). They will all be compared against monoisotopic peptide masses, but the mapping results may not be meaningful in the end." << endl;
          return false;
        }
        if (reported_mz == "average")
        {
          use_avg_mass = true;
        }
        else if (reported_mz == "maximum")
        {
          OPENMS_LOG_WARN << "For features, m/z values from the highest mass traces are reported. This type of m/z value is not available for peptides, so the comparison has to be done using average peptide masses." << endl;
          use_avg_mass = true;
        }
        before = reported_mz;
      }
    }
    return use_avg_mass;
  }

} // namespace OpenMS

