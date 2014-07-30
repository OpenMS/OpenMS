// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#ifndef OPENMS_ANALYSIS_ID_IDMAPPER_H
#define OPENMS_ANALYSIS_ID_IDMAPPER_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <limits>

namespace OpenMS
{
  /**
    @brief Annotates an MSExperiment, FeatureMap or ConsensusMap with peptide identifications

    ProteinIdentifications are assigned to the whole map.

    The retention time and mass-to-charge ratio of the PeptideIdentification as the respective "MZ" and "RT" members.

    m/z-matching on peptide side can be done either with the precursor m/z value of the peptide identification or with the theoretical masses of the peptide hits (see "mz_reference" parameter).

    See the documentation of the individual @p annotate methods for more in-depth information.

    @htmlinclude OpenMS_IDMapper.parameters

  */
  class OPENMS_DLLAPI IDMapper :
    public DefaultParamHandler
  {
public:
    enum Measure {MEASURE_PPM = 0, MEASURE_DA};

    /// Default constructor
    IDMapper();

    /// Copy C'Tor
    IDMapper(const IDMapper & cp);

    /// Assignment
    IDMapper & operator=(const IDMapper & rhs);

    /**
      @brief Mapping method for peak maps

      The identifications stored in a PeptideIdentification instance can be added to the
      corresponding spectrum.
      Note that a PeptideIdentication is added to ALL spectra which are within the allowed RT and MZ boundaries.

      @param map MSExperiment to receive the identifications
      @param ids PeptideIdentification for the MSExperiment
      @param protein_ids ProteinIdentification for the MSExperiment
      @param clear_ids Reset peptide and protein identifications of each scan before annotating
      @param mapMS1 Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)

      @exception Exception::MissingInformation is thrown if entries of @p ids do not contain 'MZ' and 'RT' information.
*/
    template <typename PeakType>
    void annotate(MSExperiment<PeakType>& map, const std::vector<PeptideIdentification>& peptide_ids, const std::vector<ProteinIdentification>& protein_ids, const bool clear_ids = false, const bool mapMS1 = false)
    {
      checkHits_(peptide_ids);

      if (clear_ids)
      { // start with empty IDs
        std::vector<PeptideIdentification> empty_ids;
        for (typename MSExperiment<PeakType>::iterator it = map.begin(); it != map.end(); ++it)
        {
          it->setPeptideIdentifications(empty_ids);
        }
        std::vector<ProteinIdentification> empty_prot_ids;
        map.setProteinIdentifications(empty_prot_ids);
      }

      if (peptide_ids.empty()) return;

      // append protein identifications
      map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

      // store mapping of scan RT to index
      std::multimap<double, Size> experiment_precursors;
      for (Size i = 0; i < map.size(); i++)
      {
        experiment_precursors.insert(std::make_pair(map[i].getRT(), i));
      }

      // store mapping of identification RT to index
      std::multimap<double, Size> identifications_precursors;
      for (Size i = 0; i < peptide_ids.size(); i++)
      {
        identifications_precursors.insert(std::make_pair(peptide_ids[i].getRT(), i));
      }

      // calculate the actual mapping
      std::multimap<double, Size>::const_iterator experiment_iterator = experiment_precursors.begin();
      std::multimap<double, Size>::const_iterator identifications_iterator = identifications_precursors.begin();
      Size matches(0);
      while (experiment_iterator != experiment_precursors.end())
      {
        // maybe we hit end() of IDs during the last scan .. go back to a real value
        if (identifications_iterator == identifications_precursors.end()) --identifications_iterator;

        // testing whether the retention times are within the precision threshold
        while (identifications_iterator != identifications_precursors.begin() &&
               (experiment_iterator->first - identifications_iterator->first) < rt_tolerance_)
        {  // go to left border of RT interval
          --identifications_iterator;
        }
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
              (identifications_iterator->first - experiment_iterator->first) < rt_tolerance_)
        {
          // testing whether the m/z fits
          if (!map[experiment_iterator->second].getPrecursors().empty() || mapMS1)
          {
            if (mapMS1 || (fabs(peptide_ids[identifications_iterator->second].getMZ() - map[experiment_iterator->second].getPrecursors()[0].getMZ()) < mz_tolerance_))
            {
              if (!(peptide_ids[identifications_iterator->second].empty()))
              {
                map[experiment_iterator->second].getPeptideIdentifications().push_back(peptide_ids[identifications_iterator->second]);
                ++matches;
              }
            }
          }
          ++identifications_iterator;
        }
        ++experiment_iterator;
      }

      // some statistics output
      LOG_INFO << "Unassigned peptides: " << peptide_ids.size() - matches << "\n"
               << "Peptides assigned to a precursor: " << matches << std::endl;

    }

    /**
      @brief Mapping method for peak maps

      Add peptide identifications stored in a feature map to their
      corresponding spectrum.

      This function converts the feature map to a vector of peptide identifications (all peptide IDs from each feature are taken)
      and calls the respective annotate() function.
      RT and m/z are taken from the peptides, or (if missing) from the feature itself.

      @param map MSExperiment to receive the identifications
      @param fmap FeatureMap with PeptideIdentifications for the MSExperiment
      @param clear_ids Reset peptide and protein identifications of each scan before annotating
      @param mapMS1 attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)

*/
    template <typename PeakType>
    void annotate(MSExperiment<PeakType>& map, FeatureMap<> fmap, const bool clear_ids = false, const bool mapMS1 = false)
    {
      const std::vector<ProteinIdentification>& protein_ids = fmap.getProteinIdentifications();
      std::vector<PeptideIdentification> peptide_ids;

      for (FeatureMap<>::const_iterator it = fmap.begin(); it != fmap.end(); ++it)
      {
        const std::vector<PeptideIdentification>& pi = it->getPeptideIdentifications();
        for (std::vector<PeptideIdentification>::const_iterator itp = pi.begin(); itp != pi.end(); ++itp)
        {
          peptide_ids.push_back(*itp);
          // if pepID has no m/z or RT, use the values of the feature
          if (!itp->hasMZ()) peptide_ids.back().setMZ(it->getMZ());
          if (!itp->hasRT()) peptide_ids.back().setRT(it->getRT());
        }

      }
      annotate(map, peptide_ids, protein_ids, clear_ids, mapMS1);
    }


    /**
      @brief Mapping method for feature maps

      If @em all features have at least one convex hull, peptide positions are matched against the bounding boxes of the convex hulls by default. If not, the positions of the feature centroids are used. The respective coordinates of the centroids are also used for matching (in place of the corresponding ranges from the bounding boxes) if @p use_centroid_rt or @p use_centroid_mz are true.

      In any case, tolerance in RT and m/z dimension is applied according to the global parameters @p rt_tolerance and @p mz_tolerance. Tolerance is understood as "plus or minus x", so the matching range is actually increased by twice the tolerance value.

      If several features (incl. tolerance) overlap the position of a peptide identification, the identification is annotated to all of them.

      @param map FeatureMap to receive the identifications
      @param ids PeptideIdentification for the ConsensusFeatures
      @param protein_ids ProteinIdentification for the ConsensusMap
      @param use_centroid_rt Whether to use the RT value of feature centroids even if convex hulls are present
      @param use_centroid_mz Whether to use the m/z value of feature centroids even if convex hulls are present

      @exception Exception::MissingInformation is thrown if entries of @p ids do not contain 'MZ' and 'RT' information.

    */
    template <typename FeatureType>
    void annotate(FeatureMap<FeatureType> & map, const std::vector<PeptideIdentification> & ids, const std::vector<ProteinIdentification> & protein_ids, bool use_centroid_rt = false, bool use_centroid_mz = false)
    {
      // std::cout << "Starting annotation..." << std::endl;
      checkHits_(ids); // check RT and m/z are present

      // append protein identifications
      map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

      // check if all features have at least one convex hull
      // if not, use the centroid and the given tolerances
      if (!(use_centroid_rt && use_centroid_mz))
      {
        for (typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it != map.end(); ++f_it)
        {
          if (f_it->getConvexHulls().empty())
          {
            use_centroid_rt = true;
            use_centroid_mz = true;
            LOG_WARN << "IDMapper warning: at least one feature has no convex hull - using centroid coordinates for matching" << std::endl;
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
      std::vector<DBoundingBox<2> > boxes;
      double min_rt = std::numeric_limits<double>::max();
      double max_rt = -std::numeric_limits<double>::max();
      // std::cout << "Precomputing bounding boxes..." << std::endl;
      boxes.reserve(map.size());
      for (typename FeatureMap<FeatureType>::Iterator f_it = map.begin();
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

        min_rt = std::min(min_rt, box.minPosition().getX());
        max_rt = std::max(max_rt, box.maxPosition().getX());
      }

      // hash bounding boxes of features by RT:
      // RT range is partitioned into slices (bins) of 1 second; every feature
      // that overlaps a certain slice is hashed into the corresponding bin
      std::vector<std::vector<SignedSize> > hash_table;
      // make sure the RT hash table has indices >= 0 and doesn't waste space
      // in the beginning:
      SignedSize offset(0);

      if (map.size() > 0)
      {
        // std::cout << "Setting up hash table..." << std::endl;
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
        LOG_WARN << "IDMapper received an empty FeatureMap! All peptides are mapped as 'unassigned'!" << std::endl;
      }

      // for statistics:
      Size matches_none = 0, matches_single = 0, matches_multi = 0;

      // std::cout << "Finding matches..." << std::endl;
      // iterate over peptide IDs:
      for (std::vector<PeptideIdentification>::const_iterator id_it =
             ids.begin(); id_it != ids.end(); ++id_it)
      {
        // std::cout << "Peptide ID: " << id_it - ids.begin() << std::endl;

        if (id_it->getHits().empty()) continue;

        DoubleList mz_values;
        double rt_value;
        IntList charges;
        getIDDetails_(*id_it, rt_value, mz_values, charges, use_avg_mass);

        if ((rt_value < min_rt) || (rt_value > max_rt))             // RT out of bounds
        {
          map.getUnassignedPeptideIdentifications().push_back(*id_it);
          ++matches_none;
          continue;
        }

        // iterate over candidate features:
        Size index = SignedSize(floor(rt_value)) - offset;
        Size matching_features = 0;
        for (std::vector<SignedSize>::iterator hash_it =
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
              for (std::vector<ConvexHull2D>::iterator ch_it =
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
                if (box.encloses(id_pos))                     // success!
                {
                  feat.getPeptideIdentifications().push_back(*id_it);
                  ++matching_features;
                  found_match = true;
                  break;                       // "ch_it" loop
                }
              }
              if (found_match) break;                   // "mz_it" loop
            }
          }
        }
        if (matching_features == 0)
        {
          map.getUnassignedPeptideIdentifications().push_back(*id_it);
          ++matches_none;
        }
        else if (matching_features == 1) ++matches_single;
        else ++matches_multi;
      }

      // some statistics output
      LOG_INFO << "Unassigned peptides: " << matches_none << "\n"
               << "Peptides assigned to exactly one feature: " << matches_single << "\n"
               << "Peptides assigned to multiple features: " << matches_multi << "\n"
               << map.getAnnotationStatistics()
               << std::endl;

    }

    /**
        @brief Mapping method for consensus maps

      If several consensus features lie inside the allowed deviation, the peptide identifications
      are mapped to all the consensus features.

      @param map ConsensusMap to receive the identifications
      @param ids PeptideIdentification for the ConsensusFeatures
      @param protein_ids ProteinIdentification for the ConsensusMap
      @param measure_from_subelements Do distance estimate from FeatureHandles instead of Centroid

        @exception Exception::MissingInformation is thrown if the MetaInfoInterface of @p ids does not contain 'MZ' and 'RT'
    */
    void annotate(ConsensusMap & map, const std::vector<PeptideIdentification> & ids, const std::vector<ProteinIdentification> & protein_ids, bool measure_from_subelements = false);

protected:
    void updateMembers_();

    /// Allowed RT deviation
    double rt_tolerance_;
    /// Allowed m/z deviation
    double mz_tolerance_;
    /// Measure used for m/z
    Measure measure_;
    /// Ignore charge states during matching?
    bool ignore_charge_;

    /// compute absolute Da tolerance, for a given m/z,
    /// when @p measure is MEASURE_DA, the value is unchanged,
    /// for MEASURE_PPM it is computed according to currently allowed ppm tolerance
    double getAbsoluteMZTolerance_(const double mz) const;

    /// check if distance constraint is fulfilled (using @p rt_tolerance_, @p mz_tolerance_ and @p measure_)
    bool isMatch_(const double rt_distance, const double mz_theoretical, const double mz_observed) const;

    /// helper function that checks if all peptide hits are annotated with RT and MZ meta values
    void checkHits_(const std::vector<PeptideIdentification> & ids) const;

    /// get RT, m/z and charge value(s) of a PeptideIdentification
    /// - multiple m/z values are returned if "mz_reference" is set to "peptide" (one for each PeptideHit)
    /// - one m/z value is returned if "mz_reference" is set to "precursor"
    void getIDDetails_(const PeptideIdentification & id, double & rt_pep, DoubleList & mz_values, IntList & charges, bool use_avg_mass = false) const;

    /// increase a bounding box by the given RT and m/z tolerances
    void increaseBoundingBox_(DBoundingBox<2> & box);

    /// try to determine the type of m/z value reported for features, return
    /// whether average peptide masses should be used for matching
    bool checkMassType_(const std::vector<DataProcessing> & processing) const;

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDMAPPER_H
