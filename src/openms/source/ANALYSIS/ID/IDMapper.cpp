// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

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
    defaults_.setValidStrings("mz_measure", ListUtils::create<String>("ppm,Da"));
    defaults_.setValue("mz_reference", "precursor", "source of m/z values for peptide identifications");
    defaults_.setValidStrings("mz_reference", ListUtils::create<String>("precursor,peptide"));

    defaults_.setValue("ignore_charge", "false", "For feature/consensus maps: Assign an ID independently of whether its charge state matches that of the (consensus) feature.");
    defaults_.setValidStrings("ignore_charge", ListUtils::create<String>("true,false"));

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

    if (clear_ids)
    { // start with empty IDs
      vector<PeptideIdentification> empty_ids;
      for (PeakMap::iterator it = map.begin(); it != map.end(); ++it)
      {
        it->setPeptideIdentifications(empty_ids);
      }
      vector<ProteinIdentification> empty_prot_ids;
      map.setProteinIdentifications(empty_prot_ids);
    }

    if (peptide_ids.empty()) return;

    // append protein identifications
    map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

    // store mapping of scan RT to index
    multimap<double, Size> experiment_precursors;
    for (Size i = 0; i < map.size(); i++)
    {
      experiment_precursors.insert(make_pair(map[i].getRT(), i));
    }

    // store mapping of identification RT to index (ignore empty hits)
    multimap<double, Size> identifications_precursors;
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      if (!peptide_ids[i].empty())
      {
        identifications_precursors.insert(make_pair(peptide_ids[i].getRT(), i));
      }
    }
    // note that mappings are sorted by key via multimap (we rely on that down below)

    // remember which peptides were mapped (for stats later)
    set<Size> peptides_mapped;

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
    std::map<Size, Size> assigned_ids;

    // keep track of assigned/unassigned precursors
    std::map<Size, Size> assigned_precursors;

    // store which peptides fit which feature (and avoid double entries)
    // consensusMap -> {peptide_index}
    vector<set<size_t>> mapping(map.size());

    DoubleList mz_values;
    double rt_pep;
    IntList charges;

    // for statistics
    Size id_matches_none(0), id_matches_single(0), id_matches_multiple(0);

    // iterate over the peptide IDs
    for (Size i = 0; i < ids.size(); ++i)
    {
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
            if (isMatch_(rt_pep - map[cm_index].getRT(), mz_pep, map[cm_index].getMZ()) && (ignore_charge_ || ListUtils::contains(current_charges, map[cm_index].getCharge())))
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
              if (isMatch_(rt_pep - it_handle->getRT(), mz_pep, it_handle->getMZ())  && (ignore_charge_ || ListUtils::contains(current_charges, it_handle->getCharge())))
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

    for (std::map<Size, Size>::const_iterator it = assigned_ids.begin(); it != assigned_ids.end(); ++it)
    {
      if (it->second == 1)
      {
        ++id_matches_single;
      }
      else if (it->second > 1)
      {
        ++id_matches_multiple;
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
          precursor_empty_id.setMetaValue("spectrum_reference",  spectra[spectrum_index].getNativeID());
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

                  // we use no undesrscore here to be compatible with linkers
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

    for (std::map<Size, Size>::const_iterator it = assigned_precursors.begin(); it != assigned_precursors.end(); ++it)
    {
      if (it->second == 1)
      {
        ++spectrum_matches_single;
      }
      else if (it->second > 1)
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

  void IDMapper::annotate(FeatureMap& map, const vector<PeptideIdentification>& ids, const vector<ProteinIdentification>& protein_ids,
                          bool use_centroid_rt, bool use_centroid_mz, const PeakMap& spectra)
  {
    // cout << "Starting annotation..." << endl;
    checkHits_(ids); // check RT and m/z are present

    // append protein identifications
    map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

    // check if all features have at least one convex hull
    // if not, use the centroid and the given tolerances
    if (!(use_centroid_rt && use_centroid_mz))
    {
      for (FeatureMap::Iterator f_it = map.begin(); f_it != map.end(); ++f_it)
      {
        if (f_it->getConvexHulls().empty())
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
    for (FeatureMap::Iterator f_it = map.begin();
         f_it != map.end(); ++f_it)
    {
      DBoundingBox<2> box;
      if (!(use_centroid_rt && use_centroid_mz))
      {
        box = f_it->getConvexHull().getBoundingBox();
      }
      if (use_centroid_rt)
      {
        box.setMinX(f_it->getRT());
        box.setMaxX(f_it->getRT());
      }
      if (use_centroid_mz)
      {
        box.setMinY(f_it->getMZ());
        box.setMaxY(f_it->getMZ());
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

    if (map.size() > 0)
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
    for (vector<PeptideIdentification>::const_iterator id_it =
         ids.begin(); id_it != ids.end(); ++id_it)
    {
      // cout << "Peptide ID: " << id_it - ids.begin() << endl;

      if (id_it->getHits().empty()) continue;

      DoubleList mz_values;
      double rt_value;
      IntList charges;
      getIDDetails_(*id_it, rt_value, mz_values, charges, use_avg_mass);

      if ((rt_value < min_rt) || (rt_value > max_rt)) // RT out of bounds
      {
        map.getUnassignedPeptideIdentifications().push_back(*id_it);
        ++matches_none;
        continue;
      }

      // iterate over candidate features:
      Size index = SignedSize(floor(rt_value)) - offset;
      Size matching_features = 0;
      for (vector<SignedSize>::iterator hash_it =
           hash_table[index].begin(); hash_it != hash_table[index].end();
           ++hash_it)
      {
        Feature & feat = map[*hash_it];

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
          if (boxes[*hash_it].encloses(id_pos))                 // potential match
          {
            if (use_centroid_mz)
            {
              // only one m/z value to check, which was already incorporated
              // into the overall bounding box -> success!
              feat.getPeptideIdentifications().push_back(*id_it);
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
                feat.getPeptideIdentifications().push_back(*id_it);
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
        map.getUnassignedPeptideIdentifications().push_back(*id_it);
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
    Size spectrum_matches(0);
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
          precursor_empty_id.setMetaValue("spectrum_reference",  spectra[spectrum_index].getNativeID());
        }
        precursor_empty_id.setIdentifier(empty_protein_id.getIdentifier());
        //precursor_empty_id.setCharge(z_p);

        for (vector<SignedSize>::iterator hash_it =
           hash_table[index].begin(); hash_it != hash_table[index].end();
           ++hash_it)
        {
          Feature & feat = map[*hash_it];

          // (optinally) check charge state
          if (!ignore_charge_)
          {
            if (z_p != feat.getCharge()) continue;
          }

          DPosition<2> id_pos(rt_value, mz_p);

          if (boxes[*hash_it].encloses(id_pos)) // potential match
          {
            if (use_centroid_mz)
            {
              // only one m/z value to check, which was already incorporated
              // into the overall bounding box -> success!
              feat.getPeptideIdentifications().push_back(precursor_empty_id);
              ++spectrum_matches;
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

    for (vector<PeptideHit>::const_iterator hit_it = id.getHits().begin();
         hit_it != id.getHits().end(); ++hit_it)
    {
      Int charge = hit_it->getCharge();
      charges.push_back(charge);

      if (param_.getValue("mz_reference") == "peptide") // use mass of each pepHit (assuming H+ adducts)
      {
        double mass = use_avg_mass ?
                      hit_it->getSequence().getAverageWeight(Residue::Full, charge) :
                      hit_it->getSequence().getMonoWeight(Residue::Full, charge);

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
    for (vector<DataProcessing>::const_iterator proc_it = processing.begin();
         proc_it != processing.end(); ++proc_it)
    {
      if (proc_it->getSoftware().getName() == "FeatureFinder")
      {
        String reported_mz = proc_it->
                             getMetaValue("parameter: algorithm:feature:reported_mz");
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
