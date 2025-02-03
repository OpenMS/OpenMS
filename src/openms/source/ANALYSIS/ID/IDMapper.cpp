// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/File.h>

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

    defaultsToParam_();
  }

  IDMapper::IDMapper(const IDMapper& cp) :
    DefaultParamHandler(cp),
    rt_tolerance_(cp.rt_tolerance_),
    mz_tolerance_(cp.mz_tolerance_),
    measure_(cp.measure_),
    ignore_charge_(cp.ignore_charge_)
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
    updateMembers_();

    return *this;
  }

  void IDMapper::updateMembers_()
  {
    rt_tolerance_ = param_.getValue("rt_tolerance");
    mz_tolerance_ = param_.getValue("mz_tolerance");
    measure_ = param_.getValue("mz_measure") == "ppm" ? MEASURE_PPM : MEASURE_DA;
    ignore_charge_ = param_.getValue("ignore_charge") == "true";
  }

  void IDMapper::annotate(PeakMap& map, const vector<PeptideIdentification>& peptide_ids, const vector<ProteinIdentification>& protein_ids, const bool clear_ids, const bool map_ms1)
  {
    checkHits_(peptide_ids);
    SpectrumLookup lookup;

    if (clear_ids)
    { // start with empty IDs
      for (PeakMap::iterator it = map.begin(); it != map.end(); ++it)
      {
        it->setPeptideIdentifications({});
      }
      map.setProteinIdentifications({});
    }

    if (peptide_ids.empty()) return;

    // append protein identifications
    map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

    lookup.readSpectra(map);

    // remember which peptides were mapped (for stats later)
    unordered_set<Size> peptides_mapped;
    // store mapping of identification RT to index (ignore empty hits)
    multimap<double, Size> identifications_precursors;
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      if (!peptide_ids[i].empty())
      { // mapping is done by either native id or by comparing peptide_id RT with experiment RT
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
            map[spectrum_idx].getPeptideIdentifications().push_back(peptide_ids[i]);
            peptides_mapped.insert(i);
          } 
          catch (const Exception::ElementNotFound& /*e*/) 
          { // use RT for mapping
            identifications_precursors.insert(make_pair(peptide_ids[i].getRT(), i));
          }
        }
      }
    }

    if (!identifications_precursors.empty()) 
    {
      // store mapping of scan RT to index
      multimap<double, Size> experiment_precursors;
      for (Size i = 0; i < map.size(); i++)
      {
        experiment_precursors.insert(make_pair(map[i].getRT(), i));
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
            for (const auto& precursor : map[experiment_iterator->second].getPrecursors())
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
            map[experiment_iterator->second].getPeptideIdentifications().push_back(peptide_ids[identifications_iterator->second]);
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


  void IDMapper::annotate(PeakMap& map, FeatureMap fmap, const bool clear_ids, const bool map_ms1)
  {
    const vector<ProteinIdentification>& protein_ids = fmap.getProteinIdentifications();
    vector<PeptideIdentification> peptide_ids;

    for (FeatureMap::const_iterator it = fmap.begin(); it != fmap.end(); ++it)
    {
      const vector<PeptideIdentification>& pi = it->getPeptideIdentifications();
      for (vector<PeptideIdentification>::const_iterator itp = pi.begin(); itp != pi.end(); ++itp)
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

  void IDMapper::annotate(
    ConsensusMap& map,
    const vector<PeptideIdentification>& ids,
    const vector<ProteinIdentification>& protein_ids,
    bool measure_from_subelements,
    bool annotate_ids_with_subelements,
    const PeakMap& spectra)
  {
    // validate "RT" and "MZ" metavalues exist
    checkHits_(ids);

    // append protein identifications to Map
    map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

    // keep track of assigned/unassigned peptide identifications.
    // maps Pep.Id. index to number of assignments to a feature
    std::unordered_map<Size, Size> assigned_ids;

    // keep track of assigned/unassigned precursors
    std::unordered_map<Size, Size> assigned_precursors;

    // store which peptides fit which feature (and avoid double entries)
    // consensusMap -> {peptide_index}
    vector<set<size_t>> mapping(map.size());

    DoubleList mz_values;
    double rt_pep;
    IntList charges;

    // for statistics
    Size id_matches_none(0), id_matches_single(0), id_matches_multiple(0);

    // build map from file to peptide id
    std::map<String, std::unordered_map<String, const PeptideIdentification*>> file2nativeid2pepid;
    bool has_spectrum_references{false};

    std::unordered_map<String, ConsensusFeature*> nativeid2cf;

    NATIVE_ID_TYPE native_id_type = checkTMTType(map);

    // We have TMT data: spectrum references annotated at consensus feature and in id
    // We can directly map by native id
    if ((native_id_type != NATIVE_ID_TYPE::UNKNOWN) )
    {
      bool lookForScanNrsAsIntegers = false;
      ProteinIdentification::Mapping mspath_mapping{protein_ids}; // used to retrieve spectrum file information annotated in protein ids given a peptide identification

      for (Size i = 0; i < ids.size(); ++i)
      {
        const PeptideIdentification* pid = &ids[i];
        
        if (pid->getHits().empty()) continue; // skip IDs without peptide annotations

        String spectrum_file = File::basename(mspath_mapping.getPrimaryMSRunPath(*pid));
        String spectrum_reference = pid->getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, "");
        // missing file origin is fine, but we need a spectrum_reference if we want to build the map
        if (spectrum_reference.empty()) continue;
        // TODO make a unique decision in the whole class on if to extract by scan number or full string?
        if (!lookForScanNrsAsIntegers)
        {
          // check if spectrum reference is a string that just contains a number
          try
          {
            ids[0].getSpectrumReference().toInt64();
            lookForScanNrsAsIntegers = true;
          }
          catch (...)
          {
            lookForScanNrsAsIntegers = false;
          }  
        }
    
        // TODO: check if there is already an entry
        file2nativeid2pepid[spectrum_file][spectrum_reference] = pid;
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

      for (auto& cf : map)
      {  
        const auto first_channel = *cf.getFeatures().begin();                  
        String filename = File::basename(map.getColumnHeaders()[first_channel.getMapIndex()].filename); // all channels are associated with same file in TMT/iTRAQ

        boost::regex scanregex{""};
        String cf_scan_id_key_name = (native_id_type == NATIVE_ID_TYPE::MS2IDMS3TMT) ? "id_scan_id" : "scan_id";
        String cf_scan_id = cf.getMetaValue(cf_scan_id_key_name, "");
        if (!cf_scan_id.empty()) 
        {
          // This assumes all scan_ids are of the same structure
          if (lookForScanNrsAsIntegers && scanregex.empty())
          {
            scanregex = SpectrumLookup::getRegExFromNativeID(cf_scan_id);
          }
          if (auto run_it = file2nativeid2pepid.find(filename); run_it != file2nativeid2pepid.end()) // TMT/iTRAQ run has identifications
          {
            if (auto scanid_it = run_it->second.find(cf_scan_id); scanid_it != run_it->second.end()) // TMT/iTRAQ run has scan_id with identification
            {
              cf.getPeptideIdentifications().push_back(*scanid_it->second);
              ++id_matches_single; // in TMT we only match to single consensus feature
            }
            // look for only the scan_number in case the search engine only extracted this (e.g. Sage)
            else if (lookForScanNrsAsIntegers)
            {
              auto scanid_it = run_it->second.find(SpectrumLookup::extractScanNumber(cf_scan_id, scanregex, false));
              if(scanid_it != run_it->second.end())
              {
                cf.getPeptideIdentifications().push_back(*scanid_it->second);
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
    }
    else
    { // non TMT data (e.g., label-free)
      for (Size i = 0; i < ids.size(); ++i)
      {
        // skip IDs without peptide annotations
        if (ids[i].getHits().empty()) continue;

        getIDDetails_(ids[i], rt_pep, mz_values, charges);

        bool id_mapped(false);

        // iterate over the features
        for (Size cm_index = 0; cm_index < map.size(); ++cm_index)
        {
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
                map[cm_index].getPeptideIdentifications().push_back(ids[i]);
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
                  if (mapping[cm_index].count(i) == 0)
                  {
                    // Store the map index of the peptide feature in the id the feature was mapped to.
                    PeptideIdentification id_pep = ids[i];
                    if (annotate_ids_with_subelements)
                    {
                      id_pep.setMetaValue("map_index", it_handle->getMapIndex());
                    }

                    map[cm_index].getPeptideIdentifications().push_back(id_pep);
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

        // the id has not been mapped to any consensus feature
        if (!id_mapped)
        {
          map.getUnassignedPeptideIdentifications().push_back(ids[i]);
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
    }

    vector<Size> unidentified = mapPrecursorsToIdentifications(spectra, ids).unidentified;

    if (!ids.empty() && !spectra.empty())
    {

      OPENMS_LOG_INFO << "Mapping " << ids.size() << "PeptideIdentifications to " << spectra.size() << " spectra." << endl;

      OPENMS_LOG_INFO << "Identification state of spectra: \n"
               << "Unidentified: " << unidentified.size() << "\n"
               << "Identified:   " << mapPrecursorsToIdentifications(spectra, ids).identified.size() << "\n"
               << "No precursor: " << mapPrecursorsToIdentifications(spectra, ids).no_precursors.size() << endl;
    }

    // we need a valid search run identifier so we try to:
    //   extract one from the map (either assigned or unassigned).
    //   or fall back to a new search run identifier.
    ProteinIdentification empty_protein_id;
    if (!unidentified.empty())
    {
      empty_protein_id.setDateTime(DateTime::now());
      if (!map.getProteinIdentifications().empty())
      {
        empty_protein_id.setIdentifier(map.getProteinIdentifications()[0].getIdentifier());
      }
      else if (!map.getUnassignedPeptideIdentifications().empty())
      {
        empty_protein_id.setIdentifier(map.getUnassignedPeptideIdentifications()[0].getIdentifier());
      }
      else
      {
        // No search run identifier given so we create a new one
        empty_protein_id.setIdentifier("UNKNOWN_SEARCH_RUN_IDENTIFIER");
        map.getProteinIdentifications().push_back(empty_protein_id);
      }
    }

    // for statistics:
    Size spectrum_matches_none(0), spectrum_matches_single(0), spectrum_matches_multiple(0);

    // are there any mapped but unidentified precursors?
    for (Size ui = 0; ui != unidentified.size(); ++ui)
    {
      Size spectrum_index = unidentified[ui];
      const MSSpectrum& spectrum = spectra[spectrum_index];
      const vector<Precursor>& precursors = spectrum.getPrecursors();

      bool precursor_mapped(false);

      // check if precursor has been identified
      for (Size i_p = 0; i_p < precursors.size(); ++i_p)
      {
        // check by precursor mass and spectrum RT
        double mz_p = precursors[i_p].getMZ();
        int z_p = precursors[i_p].getCharge();
        double rt_value = spectrum.getRT();

        PeptideIdentification precursor_empty_id;
        precursor_empty_id.setRT(rt_value);
        precursor_empty_id.setMZ(mz_p);
        precursor_empty_id.setMetaValue("spectrum_index", spectrum_index);
        if (!spectra[spectrum_index].getNativeID().empty())
        {
          precursor_empty_id.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE,  spectra[spectrum_index].getNativeID());
        }
        precursor_empty_id.setIdentifier(empty_protein_id.getIdentifier());

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
              map[cm_index].getPeptideIdentifications().push_back(precursor_empty_id);
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
                  // store the map index the precursor was mapped to
                  Size map_index = it_handle->getMapIndex();

                  // we use no underscore here to be compatible with linkers
                  precursor_empty_id.setMetaValue("map_index", map_index);
                }
                map[cm_index].getPeptideIdentifications().push_back(precursor_empty_id);
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
    if (!ids.empty())
    {
      OPENMS_LOG_INFO << "Unassigned peptides: " << id_matches_none << "\n"
               << "Peptides assigned to exactly one feature: " << id_matches_single << "\n"
               << "Peptides assigned to multiple features: " << id_matches_multiple << "\n";
    }

    if (!spectra.empty())
    {
      OPENMS_LOG_INFO << "Unassigned precursors without identification: " << spectrum_matches_none << "\n"
               << "Unidentified precursor assigned to exactly one feature: " << spectrum_matches_single << "\n"
               << "Unidentified precursor assigned to multiple features: " << spectrum_matches_multiple << "\n";
    }
  }

  void IDMapper::annotate(FeatureMap& map, 
    const vector<PeptideIdentification>& ids, 
    const vector<ProteinIdentification>& protein_ids,
    bool use_centroid_rt, 
    bool use_centroid_mz, 
    const PeakMap& spectra)
  {
    // cout << "Starting annotation..." << endl;
    checkHits_(ids); // check RT and m/z are present

    // append protein identifications
    map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

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
    // iterate over peptide IDs:
    for (const PeptideIdentification& id_it : ids)
    {
      // cout << "Peptide ID: " << id_it - ids.begin() << endl;

      if (id_it.getHits().empty()) continue;

      DoubleList mz_values;
      double rt_value;
      IntList charges;
      getIDDetails_(id_it, rt_value, mz_values, charges, use_avg_mass);

      if ((rt_value < min_rt) || (rt_value > max_rt)) // RT out of bounds
      {
        map.getUnassignedPeptideIdentifications().push_back(id_it);
        ++matches_none;
        continue;
      }

      // iterate over candidate features:
      Size index = SignedSize(floor(rt_value)) - offset;
      Size matching_features = 0;
      for (SignedSize& hash_it : hash_table[index])
      {
        Feature & feat = map[hash_it];

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
              feat.getPeptideIdentifications().push_back(id_it);
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
                feat.getPeptideIdentifications().push_back(id_it);
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
        map.getUnassignedPeptideIdentifications().push_back(id_it);
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

    vector<Size> unidentified = mapPrecursorsToIdentifications(spectra, ids).unidentified;

    // map all unidentified precursor to features
    Size spectrum_matches_none(0);
    Size spectrum_matches_single(0);
    Size spectrum_matches_multi(0);

    // we need a valid search run identifier so we try to:
    //   extract one from the map (either assigned or unassigned).
    //   or fall back to a new search run identifier.
    ProteinIdentification empty_protein_id;
    if (!unidentified.empty())
    {
      empty_protein_id.setDateTime(DateTime::now());
      if (!map.getProteinIdentifications().empty())
      {
        empty_protein_id.setIdentifier(map.getProteinIdentifications()[0].getIdentifier());
      }
      else if (!map.getUnassignedPeptideIdentifications().empty())
      {
        empty_protein_id.setIdentifier(map.getUnassignedPeptideIdentifications()[0].getIdentifier());
      }
      else
      {
        // add a new search identification run (mandatory)
        empty_protein_id.setIdentifier("UNKNOWN_SEARCH_RUN_IDENTIFIER");
        map.getProteinIdentifications().push_back(empty_protein_id);
      }
    }

    // are there any mapped but unidentified precursors?
    for (Size i = 0; i != unidentified.size(); ++i)
    {
      Size spectrum_index = unidentified[i];
      const MSSpectrum& spectrum = spectra[spectrum_index];
      const vector<Precursor>& precursors = spectrum.getPrecursors();

      // check if precursor has been identified
      for (Size i_p = 0; i_p < precursors.size(); ++i_p)
      {
        // check by precursor mass and spectrum RT
        double mz_p = precursors[i_p].getMZ();
        double rt_value = spectrum.getRT();
        int z_p = precursors[i_p].getCharge();

        if ((rt_value < min_rt) || (rt_value > max_rt)) // RT out of bounds
        {
          ++spectrum_matches_none;
          continue;
        }

        // iterate over candidate features:
        Size index = SignedSize(floor(rt_value)) - offset;
        Size matching_features = 0;

        PeptideIdentification precursor_empty_id;
        precursor_empty_id.setRT(rt_value);
        precursor_empty_id.setMZ(mz_p);
        precursor_empty_id.setMetaValue("spectrum_index", spectrum_index);
        if (!spectra[spectrum_index].getNativeID().empty())
        {
          precursor_empty_id.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, spectra[spectrum_index].getNativeID());
        }
        precursor_empty_id.setIdentifier(empty_protein_id.getIdentifier());
        //precursor_empty_id.setCharge(z_p);

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
              feat.getPeptideIdentifications().push_back(precursor_empty_id);
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
                feat.getPeptideIdentifications().push_back(precursor_empty_id);
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
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IDMapper::getAbsoluteTolerance_(): illegal internal state of measure_!", String(measure_));
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
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IDMapper::getAbsoluteTolerance_(): illegal internal state of measure_!", String(measure_));
  }

  void IDMapper::checkHits_(const vector<PeptideIdentification>& ids) const
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

  void IDMapper::getIDDetails_(const PeptideIdentification& id, double& rt_pep, DoubleList& mz_values, IntList& charges, bool use_avg_mass) const
  {
    mz_values.clear();
    charges.clear();

    rt_pep = id.getRT();

    // collect m/z values of pepId
    if (param_.getValue("mz_reference") == "precursor") // use precursor m/z of pepId
    {
      mz_values.push_back(id.getMZ());
    }

    for (const PeptideHit& hit_it : id.getHits())
    {
      Int charge = hit_it.getCharge();
      charges.push_back(charge);

      if (param_.getValue("mz_reference") == "peptide") // use mass of each pepHit (assuming H+ adducts)
      {
        double mass = use_avg_mass ?
                      hit_it.getSequence().getAverageWeight(Residue::Full, charge) :
                      hit_it.getSequence().getMonoWeight(Residue::Full, charge);

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
    String before;
    for (const DataProcessing& proc_it : processing)
    {
      if (proc_it.getSoftware().getName() == "FeatureFinder")
      {
        String reported_mz = proc_it.getMetaValue("parameter: algorithm:feature:reported_mz");
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
